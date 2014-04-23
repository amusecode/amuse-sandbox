#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#include "allvars.h"
#include "proto.h"


/*! \file density.c 
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and
 *  some auxiliary quantities are computed.  If the number of neighbours
 *  obtained falls outside the target range, the correct smoothing
 *  length is determined iteratively, if needed.
 */


#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif


/*! This function computes the local density for each active SPH particle,
 *  the number of neighbours in the current smoothing radius, and the
 *  divergence and curl of the velocity field.  The pressure is updated as
 *  well.  If a particle with its smoothing region is fully inside the
 *  local domain, it is not exported to the other processors. The function
 *  also detects particles that have a number of neighbours outside the
 *  allowed tolerance range. For these particles, the smoothing length is
 *  adjusted accordingly, and the density computation is executed again.
 *  Note that the smoothing length is not allowed to fall below the lower
 *  bound set by MinGasHsml.
 */
void density(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, n, ndone, npleft, maxfill, source, iter = 0;
  int level, ngrp, sendTask, recvTask, place, nexport;
  double dt_entr, tstart, tend, tstart_ngb = 0, tend_ngb = 0;
  double sumt, sumcomm, timengb, sumtimengb;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
  MPI_Status status;

#ifdef SPHS
  int k;
  double QVx[10],QVy[10],QVz[10],M[10][10];
  double tau, sourceterm, bit1, bit2, bit3;
  double csnd;
  double top, noise, sharp;
  int s;
  static double hubble_a, hubble_a2, fac_mu;
  double cosmofac;
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif

#ifdef SPHS
  if(All.ComovingIntegrationOn)
  {
	/* Factors for comoving integration of hydro */
  hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
  + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) +     All.OmegaLambda;
  hubble_a = All.Hubble * sqrt(hubble_a);
  hubble_a2 = All.Time * All.Time * hubble_a;
  fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;
  }
  else
  fac_mu = hubble_a = hubble_a2 = 1.0;
#endif

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      SphP[n].Left = SphP[n].Right = 0;

      if(P[n].Ti_endstep == All.Ti_Current)
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);



  /* we will repeat the whole thing for those particles where we didn't
   * find enough neighbours
   */
  do
    {
      i = 0;			/* begin with this index */
      ntotleft = ntot;		/* particles left for all tasks together */

      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
		ndone++;

		for(j = 0; j < NTask; j++)
		  Exportflag[j] = 0;

		density_evaluate(i, 0);

		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
			DensDataIn[nexport].Pos[0] = P[i].Pos[0];
			DensDataIn[nexport].Pos[1] = P[i].Pos[1];
			DensDataIn[nexport].Pos[2] = P[i].Pos[2];
			DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
			DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
			DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];
#ifdef RHOTHERM
			DensDataIn[nexport].Entropy = SphP[i].Entropy;
#endif
			DensDataIn[nexport].Hsml = SphP[i].Hsml;
			DensDataIn[nexport].Index = i;
			DensDataIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
	      }
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

	  for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	  tstart = second();

	  MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  tend = second();
	  timeimbalance += timediff(tstart, tend);


	  /* now do the particles that need to be exported */

	  for(level = 1; level < (1 << PTask); level++)
	    {
	      tstart = second();
	      for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		{
		  maxfill = 0;
		  for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* get the particles */
			  MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				       recvTask, TAG_DENS_A,
				       &DensDataGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				       MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);


	      tstart = second();
	      for(j = 0; j < nbuffer[ThisTask]; j++)
		density_evaluate(j, 1);
	      tend = second();
	      timecomp += timediff(tstart, tend);

	      /* do a block to explicitly measure imbalance */
	      tstart = second();
	      MPI_Barrier(MPI_COMM_WORLD);
	      tend = second();
	      timeimbalance += timediff(tstart, tend);

	      /* get the result */
	      tstart = second();
	      for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		{
		  maxfill = 0;
		  for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* send the results */
			  MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B,
				       &DensDataPartialResult[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
			  for(j = 0; j < nsend_local[recvTask]; j++)
			    {
 				source = j + noffset[recvTask];
				place = DensDataIn[source].Index;
				SphP[place].NumNgb += DensDataPartialResult[source].Ngb;
				SphP[place].Density += DensDataPartialResult[source].Density;
				SphP[place].DhsmlDensityFactor +=
DensDataPartialResult[source].DhsmlDensity;
#ifndef SPHS
				SphP[place].DivVel += DensDataPartialResult[source].Div;
				SphP[place].Rot[0] += DensDataPartialResult[source].Rot[0];
				SphP[place].Rot[1] += DensDataPartialResult[source].Rot[1];
				SphP[place].Rot[2] += DensDataPartialResult[source].Rot[2];
#endif
#ifdef RHOTHERM
				SphP[place].DenTherm += DensDataPartialResult[source].DenTherm;
#endif
#ifdef SPHS
			for(k = 0; k < 10; k++)
			{
				SphP[place].QVx[k] += DensDataPartialResult[source].QVx[k];
				SphP[place].QVy[k] += DensDataPartialResult[source].QVy[k];
				SphP[place].QVz[k] += DensDataPartialResult[source].QVz[k];
				SphP[place].M0[k] += DensDataPartialResult[source].M0[k];
			}
			for(k = 0; k < 6; k++)
				SphP[place].M1[k] += DensDataPartialResult[source].M1[k];
				SphP[place].M2[0] += DensDataPartialResult[source].M2[0];
				SphP[place].M2[1] += DensDataPartialResult[source].M2[1];
				SphP[place].M2[2] += DensDataPartialResult[source].M2[2];
				SphP[place].M3 += DensDataPartialResult[source].M3;
			for(k = 0; k < 6; k++)
				SphP[place].M4[k] += DensDataPartialResult[source].M4[k];
			for(k = 0; k < 5; k++)
				SphP[place].M5[k] += DensDataPartialResult[source].M5[k];
			for(k = 0; k < 4; k++)
				SphP[place].M6[k] += DensDataPartialResult[source].M6[k];
#endif
			    }
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);

	      level = ngrp - 1;
	    }

	  MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
	  for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}



      /* do final operations on results */
      tstart = second();
      for(i = 0, npleft = 0; i < N_gas; i++)
	{
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      {

		 dt_entr = (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) *
All.Timebase_interval;
		/* Conservative correction factor */ 
		SphP[i].DhsmlDensityFactor =
		  1 / (1 + SphP[i].Hsml * SphP[i].DhsmlDensityFactor / (NUMDIMS * SphP[i].Density));

#ifndef SPHS
		/* Calculate divergence and curl here [low order] */

		SphP[i].CurlVel = sqrt(SphP[i].Rot[0] * SphP[i].Rot[0] +
				       SphP[i].Rot[1] * SphP[i].Rot[1] +
				       SphP[i].Rot[2] * SphP[i].Rot[2]) / SphP[i].Density;

		SphP[i].DivVel /= SphP[i].Density;

#endif
		/* Calculate pressure */
#ifdef RHOTHERM
		SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].DenTherm, GAMMA);
#else

		SphP[i].Pressure =
		  (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
#endif
	      }


	      /* now check whether we had enough neighbours */

	      if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation) ||
		 (SphP[i].NumNgb > (All.DesNumNgb + All.MaxNumNgbDeviation)
		  && SphP[i].Hsml > (1.01 * All.MinGasHsml)))
		{
		  /* need to redo this particle */
		  npleft++;

		  if(SphP[i].Left > 0 && SphP[i].Right > 0)
		    if((SphP[i].Right - SphP[i].Left) < 1.0e-3 * SphP[i].Left)
		      {
			/* this one should be ok */
			npleft--;
			P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
			continue;
		      }

		  if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation))
		    SphP[i].Left = dmax(SphP[i].Hsml, SphP[i].Left);
		  else
		    {
		      if(SphP[i].Right != 0)
			{
			  if(SphP[i].Hsml < SphP[i].Right)
			    SphP[i].Right = SphP[i].Hsml;
			}
		      else
			SphP[i].Right = SphP[i].Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      printf
			("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (int) P[i].ID, SphP[i].Hsml, SphP[i].Left, SphP[i].Right,
			 (float) SphP[i].NumNgb, SphP[i].Right - SphP[i].Left, P[i].Pos[0], P[i].Pos[1],
			 P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(SphP[i].Right > 0 && SphP[i].Left > 0)
		    SphP[i].Hsml = pow(0.5 * (pow(SphP[i].Left, 3) + pow(SphP[i].Right, 3)), 1.0 / 3);
		  else
		    {
		      if(SphP[i].Right == 0 && SphP[i].Left == 0)
 			{
			/* This should be impossible! */
			printf("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n pos=(%g|%g|%g)\n",
i, ThisTask, (int) P[i].ID, SphP[i].Hsml, SphP[i].Left, SphP[i].Right,
(float) SphP[i].NumNgb, SphP[i].Right - SphP[i].Left, P[i].Pos[0], P[i].Pos[1],
P[i].Pos[2]);
			endrun(8188);
}

		      if(SphP[i].Right == 0 && SphP[i].Left > 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			      SphP[i].Hsml *=
				1 - (SphP[i].NumNgb -
				     All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
			    }
			  else
			    SphP[i].Hsml *= 1.26;
			}

		      if(SphP[i].Right > 0 && SphP[i].Left == 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			       SphP[i].Hsml *= 1 - (SphP[i].NumNgb -
All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;

                            }
			  else
			    SphP[i].Hsml /= 1.26;
			}
		    }


		  if(SphP[i].Hsml < All.MinGasHsml)
		    SphP[i].Hsml = All.MinGasHsml;
		}
	      else
		{
#ifdef SPHS
		/* Fill up Q and M for moment gradients */
		for(k = 0; k < 10; k++)
			{
			QVx[k] = SphP[i].QVx[k];
			QVy[k] = SphP[i].QVy[k];
			QVz[k] = SphP[i].QVz[k];
			M[0][k] = SphP[i].M0[k];
			}
		for(k = 0; k < 6; k++)
			M[1][k+4] = SphP[i].M1[k];
			M[2][5] = SphP[i].M2[0];
			M[2][6] = SphP[i].M2[1];
			M[2][9] = SphP[i].M2[2];
			M[3][6] = SphP[i].M3;
		for(k = 0; k < 6; k++)
			M[4][k+4] = SphP[i].M4[k];
		for(k = 0; k < 5; k++)
			M[5][k+5] = SphP[i].M5[k];
		for(k = 0; k < 4; k++)
			M[6][k+6] = SphP[i].M6[k];
		/* Sort out the Matrix symmetries here (not very neat) */
			M[1][0] = M[0][1];
			M[1][1] = M[0][4];
			M[1][2] = M[0][7];
			M[1][3] = M[0][8];
			M[2][0] = M[0][2];
			M[2][1] = M[1][2];
			M[2][2] = M[0][5];
			M[2][3] = M[0][9];
			M[2][4] = M[1][7];
			M[2][7] = M[1][5];
			M[2][8] = M[1][9];
			M[3][0] = M[0][3];
			M[3][1] = M[1][3];
			M[3][2] = M[2][3];
			M[3][3] = M[0][6];
			M[3][4] = M[1][8];
			M[3][5] = M[2][9];
			M[3][7] = M[1][9];
			M[3][8] = M[1][6];
			M[3][9] = M[2][6];
			M[4][0] = M[0][4];
			M[4][1] = M[1][4];
			M[4][2] = M[2][4];
			M[4][3] = M[3][4];
			M[5][0] = M[0][5];
			M[5][1] = M[1][5];
			M[5][2] = M[2][5];
			M[5][3] = M[3][5];
			M[5][4] = M[4][5];
			M[6][0] = M[0][6];
			M[6][1] = M[1][6];
			M[6][2] = M[2][6];
			M[6][3] = M[3][6];
			M[6][4] = M[4][6];
			M[6][5] = M[5][6];
			M[7][0] = M[0][7];
			M[7][1] = M[1][7];
			M[7][2] = M[2][7];
			M[7][3] = M[3][7];
			M[7][4] = M[4][7];
			M[7][5] = M[5][7];
			M[7][6] = M[6][7];
			M[7][7] = M[4][5];
			M[7][8] = M[4][9];
			M[7][9] = M[5][8];
			M[8][0] = M[0][8];
			M[8][1] = M[1][8];
			M[8][2] = M[2][8];
			M[8][3] = M[3][8];
			M[8][4] = M[4][8];
			M[8][5] = M[5][8];
			M[8][6] = M[6][8];
			M[8][7] = M[7][8];
			M[8][8] = M[4][6];
			M[8][9] = M[6][7];
			M[9][0] = M[0][9];
			M[9][1] = M[1][9];
			M[9][2] = M[2][9];
			M[9][3] = M[3][9];
			M[9][4] = M[4][9];
			M[9][5] = M[5][9];
			M[9][6] = M[6][9];
			M[9][7] = M[7][9];
			M[9][8] = M[8][9];
			M[9][9] = M[5][6];
		/* Now set up the matrix and find its inverse to give us the gradients */
		for(k = 0; k < 10; k++)
			gsl_vector_set (All.vecin, k, QVx[k]);
		for(k = 0; k < 10; k++)
			for(j = 0; j < 10; j++)
				gsl_matrix_set(All.matrixlu, k, j, M[k][j]);
			gsl_linalg_LU_decomp (All.matrixlu, All.perm, &s);
			gsl_linalg_LU_solve (All.matrixlu, All.perm, All.vecin, All.vecout);
		for(k = 0; k < 10; k++)
			SphP[i].QVx[k] = gsl_vector_get (All.vecout, k);
		for(k = 0; k < 10; k++)
			gsl_vector_set (All.vecin, k, QVy[k]);
		for(k = 0; k < 10; k++)
			for(j = 0; j < 10; j++)
				gsl_matrix_set(All.matrixlu, k, j, M[k][j]);
			gsl_linalg_LU_decomp (All.matrixlu, All.perm, &s);
			gsl_linalg_LU_solve (All.matrixlu, All.perm, All.vecin, All.vecout);
		for(k = 0; k < 10; k++)
			SphP[i].QVy[k] = gsl_vector_get (All.vecout, k);
		for(k = 0; k < 10; k++)
			gsl_vector_set (All.vecin, k, QVz[k]);
		for(k = 0; k < 10; k++)
			for(j = 0; j < 10; j++)
				gsl_matrix_set(All.matrixlu, k, j, M[k][j]);
			gsl_linalg_LU_decomp (All.matrixlu, All.perm, &s);
			gsl_linalg_LU_solve (All.matrixlu, All.perm, All.vecin, All.vecout);
		for(k = 0; k < 10; k++)
			SphP[i].QVz[k] = gsl_vector_get (All.vecout, k);
		/* Handle dissipation: */
		/* Velocity derivatives here */
		SphP[i].DivVel = -(SphP[i].QVx[1] + SphP[i].QVy[2] + SphP[i].QVz[3]);
		SphP[i].CurlVel = sqrt((SphP[i].QVz[2] - SphP[i].QVy[3])*(SphP[i].QVz[2] - SphP[i].QVy[3])+
(SphP[i].QVx[3] - SphP[i].QVz[1])*(SphP[i].QVx[3] - SphP[i].QVz[1])+
(SphP[i].QVy[1] - SphP[i].QVx[2])*(SphP[i].QVy[1] - SphP[i].QVx[2]));
		/* Viscosity */
		noise = 0.05;
		sharp = 1.0;
		csnd = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
		bit1 = SphP[i].QVx[4] + SphP[i].QVy[7] + SphP[i].QVz[8];
		bit2 = SphP[i].QVx[7] + SphP[i].QVy[5] + SphP[i].QVz[9];
		bit3 = SphP[i].QVx[8] + SphP[i].QVy[9] + SphP[i].QVz[6];
		top = sqrt(bit1*bit1 + bit2*bit2 + bit3*bit3);
		if(All.ComovingIntegrationOn)
			cosmofac = 3.0 * hubble_a2;
		else
			cosmofac = 0.0;
		sourceterm = top / (top + sharp * fabs(SphP[i].DivVel) / SphP[i].Hsml +
cosmofac / SphP[i].Hsml + noise * csnd /
(SphP[i].Hsml * fac_mu * SphP[i].Hsml));
		if(SphP[i].DivVel + cosmofac > 0)
			sourceterm = 0;
		/* Set up shear limiter */
		sourceterm *= (fabs(SphP[i].DivVel) + cosmofac) /
(fabs(SphP[i].DivVel) + cosmofac + SphP[i].CurlVel +
0.0001 * csnd / (SphP[i].Hsml * fac_mu));
		SphP[i].Alpvloc = All.ArtBulkViscConst * sourceterm;
		/* Handle dissipation decay terms */
		if(SphP[i].Alpv > SphP[i].Alpvloc)
			{
			if(All.ComovingIntegrationOn)
				tau = SphP[i].Hsml / (0.05 * SphP[i].MaxSignalVel) *
hubble_a2 * fac_mu;
			else
				tau = SphP[i].Hsml / (0.05 * SphP[i].MaxSignalVel);
			if(SphP[i].Alpv > 0.2)
				SphP[i].DtAlpv = (SphP[i].Alpvloc - SphP[i].Alpv) / tau;
			else
				SphP[i].DtAlpv = (SphP[i].Alpvloc - 0.2) / tau;
			}
		else
			{
			SphP[i].Alpv = SphP[i].Alpvloc;
			SphP[i].DtAlpv = 0;
			}
#endif
		P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
		}
	    }
	}
      tend = second();
      timecomp += timediff(tstart, tend);


      numlist = malloc(NTask * sizeof(int) * NTask);
      MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
      for(i = 0, ntot = 0; i < NTask; i++)
	ntot += numlist[i];
      free(numlist);

      if(ntot > 0)
	{
	  if(iter == 0)
	    tstart_ngb = second();

	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
      else
	tend_ngb = second();
    }
  while(ntot > 0);


  /* mark as active again */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


  /* collect some timing information */
  if(iter > 0)
    timengb = timediff(tstart_ngb, tend_ngb);
  else
    timengb = 0;

  MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
      All.CPU_EnsureNgb += sumtimengb / NTask;
    }
}



/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void density_evaluate(int target, int mode)
{
 int j, n, startnode, numngb, numngb_inbox;
 double h, h2, hinv, hinv3, hinv4;
 double rho, wk, dwk;
 double dx, dy, dz, r, r2, u, mass_j;
 double weighted_numngb, dhsmlrho;
 FLOAT *pos, *vel;

#ifndef SPHS
 double divv, rotv[3];
 double dvx, dvy, dvz;
 double fac;
#endif
#ifdef RHOTHERM
 double entropy, rhotherm;
#endif
#ifdef SPHS
 int k;
 double QVx[10], QVy[10], QVz[10];
 double M[10][10];
 double masswk,Vxmasswk,Vymasswk,Vzmasswk;
 double dx2,dy2,dz2,dxy,dxz,dyz;
#endif
 
 /* This for your choice of kernel */
#ifdef CT
 double alptrans = 0.33333333333;
 double normkern3, normkern4;
 double Norm;
#ifndef TWODIMS
 Norm = 2.53535;
#else
 Norm = 1.79360;
#endif
#endif

#ifdef HOCT
 double normkern3, normkern4;
 double alp = 0.75, bet = 0.5, nn = 6.0;
 double kap, normfac, Afac, Bfac, Pfac, Qfac;
 kap = 0.12942889;
 normfac = 15.652556;
 Afac = 5.68889;
 Bfac = -75.2;
 Pfac = -2.9888289;
 Qfac = 0.95236157;
#endif

#ifdef HOCTFAST
 double normkern3, normkern4;
 double alp = 0.75, bet = 0.5;
 double kap, normfac, Afac, Bfac, Pfac, Qfac;
 double stemp,atemp,btemp,stemp3,atemp3,btemp3;
 kap = 0.21411657;
 normfac = 6.5150497;
 Afac = 3.2;
 Bfac = -18.8;
 Pfac = -2.1542286;
 Qfac = 0.98101859;
#endif

#ifdef HIFAST
 double normkern3, normkern4;
 double alp = 0.75, bet = 0.5;
 double normfac, Afac, Bfac;
 double stemp,atemp,btemp,stemp3,atemp3,btemp3;
 normfac = 6.5264450;
 Afac = 3.2;
 Bfac = -18.8;
#endif

#ifdef TRUNCGAUSS
 double temp, normkern3, normkern4;
#endif

#ifdef HIGHSPLINE
 double normkern3, normkern4;
#endif

#ifdef LIQ
 double normkern3, normkern4;
 double alp = 0.3;
 double Afac=-1.45773,Bfac=3.79009,Cfac=-2.62391,Dfac=-0.291545;
 double Efac=0.583090,Ffac=0.650000,normfac=3.94760;
 double u2,u3,u4;
#endif

#ifdef WENDLAND
#ifdef WC2
 double normkern3 = 1.64129, normkern4 = 91.9120;
#endif
#ifdef WC4
 double normkern3 = 3.34225, normkern4 = 66.8451;
#endif
#ifdef WC6
 double normkern3 = 6.78895, normkern4 = 6.78895;
#endif
 double fivetim, u2, u3;
#endif

#ifdef QUARTIC
 double norm = 9.71405;
 double threetim, threetim2, threetim3;
 double fourtim, fourtim2, fourtim3;
#endif

#ifdef QUINTIC
 double norm = 17.4036;
 double fourtim, fourtim2, fourtim3;
 double fivetim, fivetim2, fivetim3;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h = SphP[target].Hsml;
#ifdef RHOTHERM
      entropy = SphP[target].Entropy;
#endif
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
#ifdef RHOTHERM
     entropy = DensDataGet[target].Entropy;
#endif
    }

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

#ifndef SPHS
  divv = rotv[0] = rotv[1] = rotv[2] = 0;
#endif

#ifdef RHOTHERM
  rhotherm = 0;
#endif

#ifdef SPHS
  memset(QVx,0,sizeof(QVx));
  memset(QVy,0,sizeof(QVy));
  memset(QVz,0,sizeof(QVz));
  memset(M,0,sizeof(M));
#endif

  rho = 0;
  weighted_numngb = 0;
  dhsmlrho = 0;

  startnode = All.MaxPart;
  numngb = 0;
  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      numngb++;

	      r = sqrt(r2);

	      u = r * hinv;


#ifdef WENDLAND
#ifdef WC2
	     fivetim = (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u);
	     wk = hinv3 * normkern3 * fivetim * (1.0 - u) * (3.0 + 18.0*u + 35.0*u*u);
	     dwk = -hinv4 * normkern4 * fivetim * u * (1.0 + 5.0*u);
#endif
#ifdef WC4
	     fivetim = (1.0 - u) * (1.0 - u) * (1.0 - u);
	     wk = hinv3 * normkern3 * fivetim * (1.0 - u) * (1.0 + 4.0*u);
	     dwk = -hinv4 * normkern4 * fivetim * u;
#endif
#ifdef WC6
	     fivetim = (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) *
(1.0 - u) * (1.0 - u);
	     u2 = u*u;
	     u3 = u2*u;
	     wk = hinv3 * normkern3 * fivetim * (1.0 - u) *
(1.0 + 8.0*u + 25.0*u2 + 32.0*u3);
	     dwk = -hinv4 * normkern4 * fivetim *
(22.0*u + 154.0*u2 + 352.0*u3);
#endif
#endif

#ifdef CS


	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		  dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		  dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		}


#endif
#ifdef QUARTIC
		if(u < 0.2)
		{
			threetim = (1.0 - u)*(1.0 - u)*(1.0 - u);
			threetim2 = (3.0/5.0 - u)*(3.0/5.0 - u)*(3.0/5.0 - u);
			threetim3 = (1.0/5.0 - u)*(1.0/5.0 - u)*(1.0/5.0 - u);
			fourtim = threetim*(1.0 - u);
			fourtim2 = threetim2*(3.0/5.0 - u);
			fourtim3 = threetim3*(1.0/5.0 - u);
			wk = hinv3 * norm * (fourtim - 5.0 * fourtim2 + 10.0 * fourtim3);
			dwk = hinv4 * norm * (-4.0 * threetim + 20.0 * threetim2 -
40.0 * threetim3);
		}
		else
		{
			if(u < 3.0/5.0)
			{
				threetim = (1.0 - u)*(1.0 - u)*(1.0 - u);
				threetim2 = (3.0/5.0 - u)*(3.0/5.0 - u)*(3.0/5.0 - u);
				fourtim = threetim*(1.0 - u);
				fourtim2 = threetim2*(3.0/5.0 - u);
				wk = hinv3 * norm * (fourtim - 5.0 * fourtim2);
				dwk = hinv4 * norm * (-4.0 * threetim + 20.0 * threetim2);
			}
			else
			{
				threetim = (1.0 - u)*(1.0 - u)*(1.0 - u);
				fourtim = threetim*(1.0 - u);
				wk = hinv3 * norm * (fourtim);
				dwk = hinv4 * norm * (-4.0 * threetim);
			}
		}
#endif
#ifdef QUINTIC
		if(u < 0.3333333333)
		{
			fourtim = (1.0 - u)*(1.0 - u)*(1.0 - u)*(1.0 - u);
			fivetim = fourtim * (1.0 - u);
			fourtim2 = (2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u);
			fivetim2 = fourtim2 * (2.0/3.0 - u);
			fourtim3 = (1.0/3.0 - u)*(1.0/3.0 - u)*(1.0/3.0 - u)*(1.0/3.0 - u);
			fivetim3 = fourtim3 * (1.0/3.0 - u);
			wk = hinv3 * norm * (fivetim - 6.0*fivetim2 + 15.0*fivetim3);
			dwk = hinv4 * norm * (-5.0*fourtim + 30.0*fourtim2 + 75.0*fourtim3);
		}
		else
		{
			if(u < 0.6666666666)
			{
			fourtim = (1.0 - u)*(1.0 - u)*(1.0 - u)*(1.0 - u);
			fivetim = fourtim * (1.0 - u);
			fourtim2 = (2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u);
			fivetim2 = fourtim2 * (2.0/3.0 - u);
			wk = hinv3 * norm * (fivetim - 6.0*fivetim2);
			dwk = hinv4 * norm * (-5.0*fourtim + 30.0*fourtim2);
			}
			else
			{
			fourtim = (1.0 - u)*(1.0 - u)*(1.0 - u)*(1.0 - u);
			fivetim = fourtim * (1.0 - u);
			wk = hinv3 * norm * fivetim;
			dwk = hinv4 * norm * -5.0*fourtim;
			}
		}
#endif

#ifdef CSCTMIX
		if(u < 0.5)
		{
			wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		}
		else
		{
			wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
			dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		}
#endif
	
#ifdef CT
		normkern3 = Norm * hinv3;
		normkern4 = Norm * hinv4;
		if(u < alptrans)
		{
			wk = normkern3 * (-2.0 * u + 1.2222222);
			dwk = normkern4 * (-2.0);
		}
		else
		{
			if(u < 0.5)
			{
				wk = normkern3 * (1.0 - 6.0 * u * u + 6.0 * u * u * u);
				dwk = normkern4 * (-12.0 * u + 18.0 * u * u);
			}
			else
			{
				wk = normkern3 * (2.0 * (1.0 - u) * (1.0 - u) * (1.0 - u));
				dwk = normkern4 * (-6.0 * (1.0 - u) * (1.0 - u));
			}
		}
#endif
#ifdef HOCT
		normkern3 = normfac * hinv3;
		normkern4 = normfac * hinv4;
		if(u < kap)
		{
			wk = normkern3 * (Pfac * u + Qfac);
			dwk = normkern4 * Pfac;
		}
		else
		{
			if(u < bet)
			{
			wk = normkern3 * (pow((1.0 - u),nn) + Afac * pow((alp-u),nn) + Bfac * pow((bet-u),nn));
			dwk = -normkern4 * (nn * pow((1.0 - u),nn-1.0) + nn * Afac * pow((alp-u),nn-1.0) +
nn * Bfac * pow((bet-u),nn-1.0));
			}
			else
			{
				if(u < alp)
				{
				wk = normkern3 * (pow((1.0 - u),nn) + Afac * pow((alp-u),nn));
				dwk = -normkern4 * (nn * pow((1.0 - u),nn-1.0) + nn * Afac * pow((alp-u),nn-1.0));
				}
				else
				{
				wk = normkern3 * pow((1.0 - u),nn);
				dwk = -normkern4 * nn * pow((1.0 - u),nn-1.0);
				}
			}
		}
#endif

#ifdef HOCTFAST
		normkern3 = normfac * hinv3;
		normkern4 = normfac * hinv4;
		if(u < kap)
		{
			wk = normkern3 * (Pfac * u + Qfac);
			dwk = normkern4 * Pfac;
		}
		else if(u < bet)
		{
			stemp = 1.0 - u;
			atemp = alp - u;
			btemp = bet - u;
			stemp3 = stemp*stemp*stemp;
			atemp3 = Afac*atemp*atemp*atemp;
			btemp3 = Bfac*btemp*btemp*btemp;
			wk = normkern3 * (stemp3*stemp + atemp3*atemp + btemp3*btemp);
			dwk = -normkern4 * 4.0 * (stemp3 + atemp3 + btemp3);
		}
		else if(u < alp)
		{
			stemp = 1.0 - u;
			atemp = alp - u;
			stemp3 = stemp*stemp*stemp;
			atemp3 = Afac*atemp*atemp*atemp;
			wk = normkern3 * (stemp3*stemp + atemp3*atemp);
			dwk = -normkern4 * 4.0 * (stemp3 + atemp3);
		}
		else
		{
			stemp = 1.0 - u;
			stemp3 = stemp*stemp*stemp;
			wk = normkern3 * stemp3*stemp;
			dwk = -normkern4 * 4.0 * stemp3;
		}
#endif

#ifdef LIQ
		normkern3 = normfac * hinv3;
		normkern4 = normfac * hinv4;
		if(u < alp)
		{
			wk = normkern3 * (Ffac - u);
			dwk = -normkern4;
		}
		else
		{
			u2 = u*u;
			u3 = u2*u;
			u4 = u3*u;
			wk = normkern3 * (Afac*u4 + Bfac*u3 + Cfac*u2 + Dfac*u + Efac);
			dwk = normkern4 * (4.0*Afac*u3 + 3.0*Bfac*u2 + 2.0*Cfac*u + Dfac);
		}
#endif

#ifdef HIFAST
		normkern3 = normfac * hinv3;
		normkern4 = normfac * hinv4;
		if(u < bet)
		{
			stemp = 1.0 - u;
			atemp = alp - u;
			btemp = bet - u;
			stemp3 = stemp*stemp*stemp;
			atemp3 = Afac*atemp*atemp*atemp;
			btemp3 = Bfac*btemp*btemp*btemp;
			wk = normkern3 * (stemp3*stemp + atemp3*atemp + btemp3*btemp);
			dwk = -normkern4 * 4.0 * (stemp3 + atemp3 + btemp3);
		}
		else if(u < alp)
		{
			stemp = 1.0 - u;
			atemp = alp - u;
			stemp3 = stemp*stemp*stemp;
			atemp3 = Afac*atemp*atemp*atemp;
			wk = normkern3 * (stemp3*stemp + atemp3*atemp);
			dwk = -normkern4 * 4.0 * (stemp3 + atemp3);
		}
		else
		{
			stemp = 1.0 - u;
			stemp3 = stemp*stemp*stemp;
			wk = normkern3 * stemp3*stemp;
			dwk = -normkern4 * 4.0 * stemp3;
		}
#endif

#ifdef TRUNCGAUSS
		normkern3 = All.NormKern * hinv3;
		normkern4 = All.NormKernD * hinv4;
		if(u < 1.0)
		{
			temp = exp(-(u*u/(All.Lf*All.Lf)));
			wk = normkern3 * (temp - All.StopKern);
			/* dwk = normkern4 * (-2.0 * u / All.Lf / All.Lf * temp); */
#ifdef GFORMONE
			dwk = normkern4 * (-2.0 * u / All.Lf / All.Lf * temp + 2.0 / All.Lf / All.Lf * All.StopKern);
#else
			dwk = normkern4 * (-2.0 * u / All.Lf / All.Lf * (temp - All.StopKern));
#endif
		}
		else
		{
			wk = 0;
			dwk = 0;
		}
#endif

#ifdef HIGHSPLINE
		normkern3 = All.Normspl * hinv3;
		normkern4 = All.Normspl * hinv4;
		if(u < All.betspl)
		{
			wk = normkern3 * (pow((All.rspl - u),All.nsspl) + All.Aspl*pow((All.alpspl-u),All.nsspl)+
All.Bspl*pow((All.betspl-u),All.nsspl));
			dwk = normkern4 * (-All.nsspl*pow((All.rspl - u),(All.nsspl-1)) -
All.Aspl*All.nsspl*pow((All.alpspl-u),(All.nsspl-1))-
All.Bspl*All.nsspl*pow((All.betspl-u),(All.nsspl-1)));
		}
		else
		{
			if(u < All.alpspl)
			{
				wk = normkern3 * (pow((All.rspl - u),All.nsspl) + All.Aspl*pow((All.alpspl-u),All.nsspl));
				dwk = normkern4 * (-All.nsspl*pow((All.rspl - u),(All.nsspl-1)) -
All.Aspl*All.nsspl*pow((All.alpspl-u),(All.nsspl-1)));
			}
			else
			{
				if(u < All.rspl)
				{
					wk = normkern3 * pow((All.rspl - u),All.nsspl);
					dwk = normkern4 * (-All.nsspl*pow((All.rspl - u),(All.nsspl-1)));
				}
				else
				{
					wk = 0;
					dwk = 0;
				}
			}
		}
#endif


	      mass_j = P[j].Mass;

	      dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);

	      weighted_numngb += NORM_COEFF * wk / hinv3;
		/* Calculate SPH density */
	      rho += mass_j * wk;
#ifdef RHOTHERM
	      rhotherm += pow(SphP[j].Entropy/entropy,ONE_OVER_GAMMA) * mass_j * wk;
#endif
#ifndef SPHS

	      if(r > 0)
		{
		  fac = mass_j * dwk / r;

		/* Calculate velocity divergence and curl [low order] */ 
		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];

		  divv -= fac * (dx * dvx + dy * dvy + dz * dvz);

		  rotv[0] += fac * (dz * dvy - dy * dvz);
		  rotv[1] += fac * (dx * dvz - dz * dvx);
		  rotv[2] += fac * (dy * dvx - dx * dvy);
		}

#endif
#ifdef SPHS
		/* High order Q and M matrices here for 2nd order gradients */
		masswk = P[j].Mass * wk;
		Vxmasswk = masswk * SphP[j].VelPred[0];
		Vymasswk = masswk * SphP[j].VelPred[1];
		Vzmasswk = masswk * SphP[j].VelPred[2];
		dx2 = dx * dx;
		dy2 = dy * dy;
		dz2 = dz * dz;
		dxy = dx * dy;
		dxz = dx * dz;
		dyz = dy * dz;
		QVx[0] += Vxmasswk;
		QVx[1] += Vxmasswk * dx;
		QVx[2] += Vxmasswk * dy;
		QVx[3] += Vxmasswk * dz;
		QVx[4] += Vxmasswk * dx2;
		QVx[5] += Vxmasswk * dy2;
		QVx[6] += Vxmasswk * dz2;
		QVx[7] += Vxmasswk * dxy;
		QVx[8] += Vxmasswk * dxz;
		QVx[9] += Vxmasswk * dyz;
		QVy[0] += Vymasswk;
		QVy[1] += Vymasswk * dx;
		QVy[2] += Vymasswk * dy;
		QVy[3] += Vymasswk * dz;
		QVy[4] += Vymasswk * dx2;
		QVy[5] += Vymasswk * dy2;
		QVy[6] += Vymasswk * dz2;
		QVy[7] += Vymasswk * dxy;
		QVy[8] += Vymasswk * dxz;
		QVy[9] += Vymasswk * dyz;
		QVz[0] += Vzmasswk;
		QVz[1] += Vzmasswk * dx;
		QVz[2] += Vzmasswk * dy;
		QVz[3] += Vzmasswk * dz;
		QVz[4] += Vzmasswk * dx2;
		QVz[5] += Vzmasswk * dy2;
		QVz[6] += Vzmasswk * dz2;
		QVz[7] += Vzmasswk * dxy;
		QVz[8] += Vzmasswk * dxz;
		QVz[9] += Vzmasswk * dyz;
		/* There's a lot of symmetry in this 10x10 matrix so
we don't need to calculate every term */
		M[0][0] += masswk;
		M[0][1] += masswk * dx;
		M[0][2] += masswk * dy;
		M[0][3] += masswk * dz;
		M[0][4] += masswk * dx2;
		M[0][5] += masswk * dy2;
		M[0][6] += masswk * dz2;
		M[0][7] += masswk * dxy;
		M[0][8] += masswk * dxz;
		M[0][9] += masswk * dyz;
		M[1][4] += masswk * dx2 * dx;
		M[1][5] += masswk * dy2 * dx;
		M[1][6] += masswk * dz2 * dx;
		M[1][7] += masswk * dxy * dx;
		M[1][8] += masswk * dxz * dx;
		M[1][9] += masswk * dyz * dx;
		M[2][5] += masswk * dy2 * dy;
		M[2][6] += masswk * dz2 * dy;
		M[2][9] += masswk * dyz * dy;
		M[3][6] += masswk * dz2 * dz;
		M[4][4] += masswk * dx2 * dx2;
		M[4][5] += masswk * dy2 * dx2;
		M[4][6] += masswk * dz2 * dx2;
		M[4][7] += masswk * dxy * dx2;
		M[4][8] += masswk * dxz * dx2;
		M[4][9] += masswk * dyz * dx2;
		M[5][5] += masswk * dy2 * dy2;
		M[5][6] += masswk * dz2 * dy2;
		M[5][7] += masswk * dxy * dy2;
		M[5][8] += masswk * dxz * dy2;
		M[5][9] += masswk * dyz * dy2;
		M[6][6] += masswk * dz2 * dz2;
		M[6][7] += masswk * dxy * dz2;
		M[6][8] += masswk * dxz * dz2;
		M[6][9] += masswk * dyz * dz2;
#endif

	    }
	}
    }
  while(startnode >= 0);

  if(mode == 0)
    {
 	SphP[target].NumNgb = weighted_numngb;
	SphP[target].Density = rho;
	SphP[target].DhsmlDensityFactor = dhsmlrho;
#ifndef SPHS
	SphP[target].DivVel = divv;
	SphP[target].Rot[0] = rotv[0];
	SphP[target].Rot[1] = rotv[1];
	SphP[target].Rot[2] = rotv[2];
#endif
#ifdef RHOTHERM
	SphP[target].DenTherm = rhotherm;
#endif
#ifdef SPHS
	for(k = 0; k < 10; k++)
	{
		SphP[target].QVx[k] = QVx[k];
		SphP[target].QVy[k] = QVy[k];
		SphP[target].QVz[k] = QVz[k];
		SphP[target].M0[k] = M[0][k];
	}
	for(k = 0; k < 6; k++)
		SphP[target].M1[k] = M[1][k+4];
	SphP[target].M2[0] = M[2][5];
	SphP[target].M2[1] = M[2][6];
	SphP[target].M2[2] = M[2][9];
	SphP[target].M3 = M[3][6];
	for(k = 0; k < 6; k++)
		SphP[target].M4[k] = M[4][k+4];
	for(k = 0; k < 5; k++)
		SphP[target].M5[k] = M[5][k+5];
	for(k = 0; k < 4; k++)
		SphP[target].M6[k] = M[6][k+6];
#endif
    }
  else
    {
	DensDataResult[target].Density = rho;
	DensDataResult[target].Ngb = weighted_numngb;
	DensDataResult[target].DhsmlDensity = dhsmlrho;
#ifndef SPHS
	DensDataResult[target].Div = divv;
	DensDataResult[target].Rot[0] = rotv[0];
	DensDataResult[target].Rot[1] = rotv[1];
	DensDataResult[target].Rot[2] = rotv[2];
#endif

#ifdef RHOTHERM
	DensDataResult[target].DenTherm = rhotherm;
#endif

#ifdef SPHS
	for(k = 0; k < 10; k++)
	{
		DensDataResult[target].QVx[k] = QVx[k];
		DensDataResult[target].QVy[k] = QVy[k];
		DensDataResult[target].QVz[k] = QVz[k];
		DensDataResult[target].M0[k] = M[0][k];
	}
	for(k = 0; k < 6; k++)
		DensDataResult[target].M1[k] = M[1][k+4];
	DensDataResult[target].M2[0] = M[2][5];
	DensDataResult[target].M2[1] = M[2][6];
	DensDataResult[target].M2[2] = M[2][9];
	DensDataResult[target].M3 = M[3][6];
	for(k = 0; k < 6; k++)
		DensDataResult[target].M4[k] = M[4][k+4];
	for(k = 0; k < 5; k++)
		DensDataResult[target].M5[k] = M[5][k+5];
	for(k = 0; k < 4; k++)
		DensDataResult[target].M6[k] = M[6][k+6];
#endif
    }
}




/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int dens_compare_key(const void *a, const void *b)
{
  if(((struct densdata_in *) a)->Task < (((struct densdata_in *) b)->Task))
    return -1;

  if(((struct densdata_in *) a)->Task > (((struct densdata_in *) b)->Task))
    return +1;

  return 0;
}
