#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#include "allvars.h"
#include "proto.h"

/*! \file hydra.c
 *  \brief Computation of SPH forces and rate of entropy generation
 *
 *  This file contains the "second SPH loop", where the SPH forces are
 *  computed, and where the rate of change of entropy due to the shock heating
 *  (via artificial viscosity) is computed.
 */

static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;

#ifdef SPHS
static double asqmu;
#endif

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


/*! This function is the driver routine for the calculation of hydrodynamical
 *  force and rate of change of entropy due to shock heating for all active
 *  particles .
 */
void hydro_force(void)
{
  long long ntot, ntotleft;
  int i, j, k, n, ngrp, maxfill, source, ndone;
  int *nbuffer, *noffset, *nsend_local, *nsend, *numlist, *ndonelist;
  int level, sendTask, recvTask, nexport, place;
  double soundspeed_i;
  double tstart, tend, sumt, sumcomm;
  double timecomp = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  MPI_Status status;

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

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda;

      hubble_a = All.Hubble * sqrt(hubble_a);
      hubble_a2 = All.Time * All.Time * hubble_a;

      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;

      fac_egy = pow(All.Time, 3 * (GAMMA - 1));

      fac_vsic_fix = hubble_a * pow(All.Time, 3 * GAMMA_MINUS1);

      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;

#ifdef SPHS
      asqmu = All.Time * All.Time * fac_mu;
#endif
    }
  else
    {
      hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = fac_egy = 1.0;
#ifdef SPHS
      asqmu = 1.0; 
#endif
    }

  /* `NumSphUpdate' gives the number of particles on this processor that want a force update */
  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      if(P[n].Ti_endstep == All.Ti_Current)
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;			/* first particle for this task */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeHydro - NTask; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;

	    hydro_evaluate(i, 0);

	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++)
		      {
			HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
			HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
		      }
		    HydroDataIn[nexport].Hsml = SphP[i].Hsml;
		    HydroDataIn[nexport].Mass = P[i].Mass;
		    HydroDataIn[nexport].DhsmlDensityFactor = SphP[i].DhsmlDensityFactor;
		    HydroDataIn[nexport].Density = SphP[i].Density;
		    HydroDataIn[nexport].Pressure = SphP[i].Pressure;
		    HydroDataIn[nexport].Timestep = P[i].Ti_endstep - P[i].Ti_begstep;
		    HydroDataIn[nexport].Entropy = SphP[i].Entropy;

#ifdef SPHS
                    HydroDataIn[nexport].Alpv = SphP[i].Alpv;
#endif

		    /* calculation of F1 */
		    soundspeed_i = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);

#ifndef SPHS
		    HydroDataIn[nexport].F1 = fabs(SphP[i].DivVel) /
		      (fabs(SphP[i].DivVel) + SphP[i].CurlVel +
		       0.0001 * soundspeed_i / SphP[i].Hsml / fac_mu);
#endif		
    
		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;
		    nexport++;
		    nsend_local[j]++;
		  }
	      }
	  }
      tend = second();
      timecomp += timediff(tstart, tend);

      qsort(HydroDataIn, nexport, sizeof(struct hydrodata_in), 
	    hydro_compare_key);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, 
		    MPI_COMM_WORLD);

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
	      if(maxfill >= All.BunchSizeHydro)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || 
		     nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * 
				   sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A,
				   &HydroDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * 
				   sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, 
				   &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  /* now do the imported particles */
	  tstart = second();
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    hydro_evaluate(j, 1);
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  /* do a block to measure imbalance */
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
	      if(maxfill >= All.BunchSizeHydro)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || 
		     nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  source = j + noffset[recvTask];
			  place = HydroDataIn[source].Index;

			  for(k = 0; k < 3; k++)
			    SphP[place].HydroAccel[k] += HydroDataPartialResult[source].Acc[k];
#ifdef ERRORS
                          for(k = 0; k < 3; k++)
			    SphP[place].Error0[k] += HydroDataPartialResult[source].Error0[k];
#endif
#ifdef SPHS
			  SphP[place].EntDiss += HydroDataPartialResult[source].EntDiss;
                          SphP[place].MassDiss += HydroDataPartialResult[source].MassDiss;
#endif
			  SphP[place].DtEntropy += HydroDataPartialResult[source].DtEntropy;

			  if(SphP[place].MaxSignalVel < HydroDataPartialResult[source].MaxSignalVel)
			    SphP[place].MaxSignalVel = HydroDataPartialResult[source].MaxSignalVel;
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

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);

  /* do final operations on results */
  tstart = second();

  for(i = 0; i < N_gas; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
	/* Multiply entropy by pre-factors and handle special cases */
	SphP[i].DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].Density, GAMMA_MINUS1));

#ifdef SPHS
        /* Handle dissipation terms */
	SphP[i].EntDiss *= asqmu / hubble_a2;
	SphP[i].MassDiss *= asqmu / hubble_a2;

	SphP[i].DtEntropy += SphP[i].EntDiss - 
	  SphP[i].MassDiss / P[i].Mass * fac_mu * fac_mu * 
	  (-0.5*(P[i].Vel[0]*P[i].Vel[0]+
		 P[i].Vel[1]*P[i].Vel[1]+
		 P[i].Vel[2]*P[i].Vel[2]) * GAMMA_MINUS1 / 
	   (pow(SphP[i].Density, GAMMA_MINUS1))) - 
	  SphP[i].MassDiss / P[i].Mass * SphP[i].Entropy;
	
	for(k = 0; k < 3; k++)
	  SphP[i].HydroAccel[k] -= hubble_a2 * fac_mu * fac_mu *
	    SphP[i].MassDiss / P[i].Mass * P[i].Vel[k];
#endif

	/* Special cases ... */ 
#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
	    SphP[i].DtEntropy = 0;
	    for(k = 0; k < 3; k++)
	      SphP[i].HydroAccel[k] = 0;
	  }
#endif
      }
  
  tend = second();
  timecomp += timediff(tstart, tend);

  /* collect some timing information */
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
    }
}


/*! This function is the 'core' of the SPH force computation. A target
 *  particle is specified which may either be local, or reside in the
 *  communication buffer.
 */
void hydro_evaluate(int target, int mode)
{
  int j, k, n, timestep, startnode, numngb;
  FLOAT *pos, *vel;
  FLOAT mass, h_i, dhsmlDensityFactor, rho, pressure;
  double acc[3], dtEntropy;
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv4;
  double soundspeed_i, soundspeed_j;
  double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
  double h_j, dwk_j, r, r2, u, hfc_visc;
  double entropy;
  double rhoj, pressurej;
  double maxSignalVel;
  double p_over_rho2_i, p_over_rho2_j;

#ifndef SPHS
  double f1, f2;
#endif
#ifdef SPHS
  double mvsig, alpv;
  double alpe_ij, alpm_ij, alpv_ij;
  double entdiss, massdiss;
  double pressurelimiter;
#endif

  /* This for your choice of kernel */
#ifdef CS
  double wk_i = 0, wk_j = 0, hinv3;
#endif
#ifdef CT
  double alptrans = 0.33333333333;
  double normkern4;
  double Norm;
#ifndef TWODIMS
  Norm = 2.53535;
#else
  Norm = 1.79360;
#endif
#endif

#ifdef CSCTMIX
  double alptrans = 0.33333333333;
  double normkern4;
  double Norm;
#ifndef TWODIMS
  Norm = 2.53535;
#else
  Norm = 1.79360;
#endif
#endif

#ifdef HOCT
  double normkern4;
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
  double normkern4;
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
  double normkern4;
  double alp = 0.75, bet = 0.5;
  double normfac, Afac, Bfac;
  double stemp,atemp,btemp,stemp3,atemp3,btemp3;
  normfac = 6.5264450;
  Afac = 3.2;
  Bfac = -18.8;
#endif

#ifdef TRUNCGAUSS
  double temp, normkern4;
#endif

#ifdef HIGHSPLINE
  double normkern4;
#endif

#ifdef LIQ
  double normkern4;
  double alp = 0.3;
  double Afac=-1.45773,Bfac=3.79009,Cfac=-2.62391,Dfac=-0.291545;
  double Efac=0.583090,Ffac=0.650000,normfac=3.94760;
  double u2,u3;
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
  double u2,u3;
#endif
  double fivetim;
#endif

#ifdef QUARTIC
  double norm = 9.71405;
  double threetim, threetim2, threetim3;
  double fourtim, fourtim2, fourtim3;
#endif

#ifdef QUINTIC
  double norm = 17.4036;
  double fourtim, fourtim2, fourtim3;
#endif

  /* This for calculating SPH errors */ 
#ifdef ERRORS
  double error0x,error0y,error0z;
  double phi, phij; 
#endif

#ifndef NOVISCOSITYLIMITER
  double dt;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = SphP[target].Hsml;
      mass = P[target].Mass;
      dhsmlDensityFactor = SphP[target].DhsmlDensityFactor;
      entropy = SphP[target].Entropy;
      timestep = P[target].Ti_endstep - P[target].Ti_begstep;
      rho = SphP[target].Density;
#ifdef SPHS
      alpv = SphP[target].Alpv;
#endif
      pressure = SphP[target].Pressure;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#ifndef SPHS
      f1 = fabs(SphP[target].DivVel) /
	(fabs(SphP[target].DivVel) + SphP[target].CurlVel +
	 0.0001 * soundspeed_i / SphP[target].Hsml / fac_mu);
#endif
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;
      mass = HydroDataGet[target].Mass;
      dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
      entropy = HydroDataGet[target].Entropy;
      timestep = HydroDataGet[target].Timestep;
      rho = HydroDataGet[target].Density;
#ifdef SPHS
      alpv = HydroDataGet[target].Alpv;
#endif
      pressure = HydroDataGet[target].Pressure;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#ifndef SPHS
      f1 = HydroDataGet[target].F1;
#endif
    }

  /* initialize variables before SPH loop is started */
  acc[0] = acc[1] = acc[2] = dtEntropy = 0;
  maxSignalVel = 0;

#ifdef ERRORS
  error0x = error0y = error0z = 0;
#endif
#ifdef SPHS
  entdiss = 0;
  massdiss = 0;
#endif

  p_over_rho2_i = pressure / (rho * rho) * dhsmlDensityFactor;
  h_i2 = h_i * h_i;

  /* Now start the actual SPH computation for this particle */
  startnode = All.MaxPart;
  do
    {
      numngb = ngb_treefind_pairs(&pos[0], h_i, &startnode);

      for(n = 0; n < numngb; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  find the closest image in the given box size  */
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
	  h_j = SphP[j].Hsml;

	  if(r2 < h_i2 || r2 < h_j * h_j)
	    {
	      r = sqrt(r2);
	      if(r > 0)
		{
		  rhoj = SphP[j].Density;
		  pressurej = SphP[j].Pressure;
		  soundspeed_j = sqrt(GAMMA * pressurej / rhoj);
		  p_over_rho2_j = pressurej / (rhoj * rhoj);

		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];

                  vdotr = dx * dvx + dy * dvy + dz * dvz;

		  if(All.ComovingIntegrationOn)
		    vdotr2 = vdotr + hubble_a2 * r2;
		  else
		    vdotr2 = vdotr;

#ifdef CS
		  if(r2 < h_i2)
		    {
		      hinv = 1.0 / h_i;
#ifndef  TWODIMS
		      hinv3 = hinv * hinv * hinv;
		      hinv4 = hinv3 * hinv;
#else
		      hinv3 = hinv * hinv;
		      hinv4 = hinv3 * hinv / boxSize_Z;
#endif

		      u = r * hinv;
		      
		      if(u < 0.5)
			{
			  wk_i = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			  dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			}
		      else
			{
			  wk_i = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
			  dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
			}
		    }
		  else
		    {
		      wk_i = 0;
		      dwk_i = 0;
		    }
		  
		  if(r2 < h_j * h_j)
		    {
		      hinv = 1.0 / h_j;
#ifndef  TWODIMS
                      hinv3 = hinv * hinv * hinv;
		      hinv4 = hinv3 * hinv;
#else
                      hinv3 = hinv * hinv;
		      hinv4 = hinv3 * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      if(u < 0.5)
			{
			  wk_j = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			  dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			}
		      else
			{
			  wk_j = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
			  dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
			}
		    }
		  else
		    {
		      wk_j = 0;
		      dwk_j = 0;
		    }
#endif

#ifdef QUARTIC
                  if(r2 < h_i2)
                    {
                      hinv = 1.0 / h_i;
                      hinv4 = hinv * hinv * hinv * hinv;
                      u = r * hinv;

                      if(u < 0.2)
			{
			  threetim = (1.0 - u)*(1.0 - u)*(1.0 - u);
			  threetim2 = (3.0/5.0 - u)*(3.0/5.0 - u)*(3.0/5.0 - u);
			  threetim3 = (1.0/5.0 - u)*(1.0/5.0 - u)*(1.0/5.0 - u);

			  fourtim = threetim*(1.0 - u);
			  fourtim2 = threetim2*(3.0/5.0 - u);
			  fourtim3 = threetim3*(1.0/5.0 - u);

			  dwk_i = hinv4 * norm * (-4.0 * threetim + 20.0 * threetim2 -
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
			      dwk_i = hinv4 * norm * (-4.0 * threetim + 20.0 * threetim2);
			    }
			  else
			    {
			      threetim = (1.0 - u)*(1.0 - u)*(1.0 - u);
			      fourtim = threetim*(1.0 - u);
			      dwk_i = hinv4 * norm * (-4.0 * threetim);
			    }
			}
		    }
		  else
		    {
		      dwk_i = 0;
		    }

		  if(r2 < h_j * h_j)
		    {
		      hinv = 1.0 / h_j;
		      hinv4 = hinv * hinv * hinv * hinv;
		      u = r * hinv;

		      if(u < 0.2)
			{
			  threetim = (1.0 - u)*(1.0 - u)*(1.0 - u);
			  threetim2 = (3.0/5.0 - u)*(3.0/5.0 - u)*(3.0/5.0 - u);
			  threetim3 = (1.0/5.0 - u)*(1.0/5.0 - u)*(1.0/5.0 - u);
			  
			  fourtim = threetim*(1.0 - u);
			  fourtim2 = threetim2*(3.0/5.0 - u);
			  fourtim3 = threetim3*(1.0/5.0 - u);

			  dwk_j = hinv4 * norm * (-4.0 * threetim + 20.0 * threetim2 -
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
			      dwk_j = hinv4 * norm * (-4.0 * threetim + 20.0 * threetim2);
			    }
			  else
			    {
			      threetim = (1.0 - u)*(1.0 - u)*(1.0 - u);
			      fourtim = threetim*(1.0 - u);
			      dwk_j = hinv4 * norm * (-4.0 * threetim);
			    }
			}
		    }
		  else
		    {
                      dwk_j = 0;
		    }
#endif

#ifdef QUINTIC
                  if(r2 < h_i2)
                    {
                      hinv = 1.0 / h_i;
                      hinv4 = hinv * hinv * hinv * hinv;
                      u = r * hinv;

		      if(u < 0.3333333333)
			{
			  fourtim = (1.0 - u)*(1.0 - u)*(1.0 - u)*(1.0 - u);
			  fourtim2 = (2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u);
			  fourtim3 = (1.0/3.0 - u)*(1.0/3.0 - u)*(1.0/3.0 - u)*(1.0/3.0 - u);
			  dwk_i = hinv4 * norm * (-5.0*fourtim + 30.0*fourtim2 + 75.0*fourtim3);
			}
		      else
			{
			  if(u < 0.6666666666)
			    {
			      fourtim = (1.0 - u)*(1.0 - u)*(1.0 - u)*(1.0 - u);
			      fourtim2 = (2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u);
			      dwk_i = hinv4 * norm * (-5.0*fourtim + 30.0*fourtim2);
			    }
			  else
			    {
			      fourtim = (1.0 - u)*(1.0 - u)*(1.0 - u)*(1.0 - u);
			      dwk_i = hinv4 * norm * -5.0*fourtim;
			    }
			}
		    }
		  else
		    {
		      dwk_i = 0;
		    }
		  
		  if(r2 < h_j * h_j)
		    {
		      hinv = 1.0 / h_j;
		      hinv4 = hinv * hinv * hinv * hinv;
		      u = r * hinv;
		      
		      if(u < 0.3333333333)
			{
			  fourtim = (1.0 - u)*(1.0 - u)*(1.0 - u)*(1.0 - u);
			  fourtim2 = (2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u);
			  fourtim3 = (1.0/3.0 - u)*(1.0/3.0 - u)*(1.0/3.0 - u)*(1.0/3.0 - u);
			  dwk_j = hinv4 * norm * (-5.0*fourtim + 30.0*fourtim2 + 75.0*fourtim3);
			}
		      else
			{
			  if(u < 0.6666666666)
			    {
			      fourtim = (1.0 - u)*(1.0 - u)*(1.0 - u)*(1.0 - u);
			      fourtim2 = (2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u)*(2.0/3.0 - u);
			      dwk_j = hinv4 * norm * (-5.0*fourtim + 30.0*fourtim2);
			    }
			  else
			    {
			      fourtim = (1.0 - u)*(1.0 - u)*(1.0 - u)*(1.0 - u);
			      dwk_j = hinv4 * norm * -5.0*fourtim;
			    }
			}
		    }
		  else
		    {
		      dwk_j = 0;
		    }
#endif
		  
#ifdef CT
		  if(r2 < h_i2)
		    {
		      
                      hinv = 1.0 / h_i;
#ifndef  TWODIMS
                      hinv4 = hinv * hinv * hinv * hinv;
#else
                      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      normkern4 = Norm * hinv4;

		      if(u < alptrans)
                        {
                          dwk_i = normkern4 * (-2.0);
                        }
		      else
			{
			  if(u < 0.5)
			    {
			      dwk_i = normkern4 * (-12.0 * u + 18.0 * u * u);
			    }
			  else
			    {
			      dwk_i = normkern4 * (-6.0 * (1.0 - u) * (1.0 - u));
			    }
			}
		    }
		  else
		    {
		      dwk_i = 0;
		    }
		  
		  if(r2 < h_j * h_j)
		    {
                      hinv = 1.0 / h_j;
#ifndef  TWODIMS
                      hinv4 = hinv * hinv * hinv * hinv;
#else
                      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      normkern4 = Norm * hinv4;
		      if(u < alptrans)
                        {
                          dwk_j = normkern4 * (-2.0);
                        }
		      else
			{
			  if(u < 0.5)
			    {
			      dwk_j = normkern4 * (-12.0 * u + 18.0 * u * u);
			    }
			  else
			    {
			      dwk_j = normkern4 * (-6.0 * (1.0 - u) * (1.0 - u));
			    }
			}
		    }
		  else
		    {
		      dwk_j = 0;
		    }
#endif

#ifdef CSCTMIX
		  if(r2 < h_i2)
		    {

                      hinv = 1.0 / h_i;
#ifndef  TWODIMS
                      hinv4 = hinv * hinv * hinv * hinv;
#else
                      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      normkern4 = Norm * hinv4;

		      if(u < alptrans)
			dwk_i = normkern4 * (-2.0);
		      else if(u < 0.5)
			dwk_i = normkern4 * (-12.0 * u + 18.0 * u * u);
		      else
			dwk_i = normkern4 * (-6.0 * (1.0 - u) * (1.0 - u));
		    }
		  else
		    dwk_i = 0;
		    
		  if(r2 < h_j * h_j)
		    {
                      hinv = 1.0 / h_j;
#ifndef  TWODIMS
                      hinv4 = hinv * hinv * hinv * hinv;
#else
                      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      normkern4 = Norm * hinv4;
		      if(u < alptrans)
			dwk_j = normkern4 * (-2.0);
		      else if(u < 0.5)
			dwk_j = normkern4 * (-12.0 * u + 18.0 * u * u);
		      else
			dwk_j = normkern4 * (-6.0 * (1.0 - u) * (1.0 - u));
		    }
		  else
		    dwk_j = 0;
#endif

#ifdef HOCT
                  if(r2 < h_i2)
                    {
                      hinv = 1.0 / h_i;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
                      normkern4 = normfac * hinv4;

		      if(u < kap)
			{
                          dwk_i = normkern4 * Pfac;
			}
		      else
			{
			  if(u < bet)
			    {
			      dwk_i = -normkern4 * (nn * pow((1.0 - u),nn-1.0) + nn * Afac * pow((alp-u),nn-1.0) + nn * Bfac * pow((bet-u),nn-1.0));
			    }
			  else
			    {
			      if(u < alp)
				{
				  dwk_i = -normkern4 * (nn * pow((1.0 - u),nn-1.0) + nn * Afac * pow((alp-u),nn-1.0));
				}
			      else
				{
				  dwk_i = -normkern4 * nn * pow((1.0 - u),nn-1.0);
				}
			    }
			}
		    }
		  else
		    {
		      dwk_i = 0;
		    }
		  
		  if(r2 < h_j * h_j)
		    {
                      hinv = 1.0 / h_j;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
                      normkern4 = normfac * hinv4;

                      if(u < kap)
                        {
                          dwk_j = normkern4 * Pfac;
                        }
                      else
                        {
                          if(u < bet)
                            {
                              dwk_j = -normkern4 * (nn * pow((1.0 - u),nn-1.0) + nn * Afac * pow((alp-u),nn-1.0) + nn * Bfac * pow((bet-u),nn-1.0));
                            }
                          else
                            {
                              if(u < alp)
                                {
                                  dwk_j = -normkern4 * (nn * pow((1.0 - u),nn-1.0) + nn * Afac * pow((alp-u),nn-1.0));
                                }
                              else
                                {
                                  dwk_j = -normkern4 * nn * pow((1.0 - u),nn-1.0);
                                }
                            }
                        }
                    }
                  else
                    {
                      dwk_j = 0;
                    }
#endif

#ifdef HOCTFAST
                  if(r2 < h_i2)
                    {
                      hinv = 1.0 / h_i;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
		      normkern4 = normfac * hinv4;
		      if(u < kap)
			{
			  dwk_i = normkern4 * Pfac;
			}
		      else if(u < bet)
			{
			  stemp = 1.0 - u;
			  atemp = alp - u;
			  btemp = bet - u;
			  stemp3 = stemp*stemp*stemp;
			  atemp3 = Afac*atemp*atemp*atemp;
			  btemp3 = Bfac*btemp*btemp*btemp;
			  dwk_i = -normkern4 * 4.0 * (stemp3 + atemp3 + btemp3); 
			}
		      else if(u < alp)
			{
			  stemp = 1.0 - u;
			  atemp = alp - u;
			  stemp3 = stemp*stemp*stemp;
			  atemp3 = Afac*atemp*atemp*atemp;
			  dwk_i = -normkern4 * 4.0 * (stemp3 + atemp3); 
			}
		      else
			{
			  stemp = 1.0 - u;
			  stemp3 = stemp*stemp*stemp;
			  dwk_i = -normkern4 * 4.0 * stemp3; 
			}
		    }
                  else
                    {
                      dwk_i = 0;
                    }

                  if(r2 < h_j * h_j)
                    {
                      hinv = 1.0 / h_j;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
                      normkern4 = normfac * hinv4;
		      
		      if(u < kap)
			{
			  dwk_j = normkern4 * Pfac;
			}
		      else if(u < bet)
			{
			  stemp = 1.0 - u;
			  atemp = alp - u;
			  btemp = bet - u;
			  stemp3 = stemp*stemp*stemp;
			  atemp3 = Afac*atemp*atemp*atemp;
			  btemp3 = Bfac*btemp*btemp*btemp;
			  dwk_j = -normkern4 * 4.0 * (stemp3 + atemp3 + btemp3); 
			}
		      else if(u < alp)
			{
			  stemp = 1.0 - u;
			  atemp = alp - u;
			  stemp3 = stemp*stemp*stemp;
			  atemp3 = Afac*atemp*atemp*atemp;
			  dwk_j = -normkern4 * 4.0 * (stemp3 + atemp3);
			}
		      else
			{
			  stemp = 1.0 - u;
			  stemp3 = stemp*stemp*stemp;
			  dwk_j = -normkern4 * 4.0 * stemp3;
			}
		    }
                  else
                    {
                      dwk_j = 0; 
                    }
#endif

#ifdef LIQ
                  if(r2 < h_i2)
                    {
                      hinv = 1.0 / h_i;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
		      normkern4 = normfac * hinv4;
		      if(u < alp)
			{
			  dwk_i = -normkern4;
			}
		      else
			{
			  u2 = u*u;
			  u3 = u2*u;
			  dwk_i = normkern4 * (4.0*Afac*u3 + 3.0*Bfac*u2 + 2.0*Cfac*u + Dfac);
			}
		    }
                  else
                    {
                      dwk_i = 0;
                    }

                  if(r2 < h_j * h_j)
                    {
                      hinv = 1.0 / h_j;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
                      normkern4 = normfac * hinv4;
		            
		      if(u < alp)
			{
			  dwk_j = -normkern4;
			}
		      else
			{
			  u2 = u*u;
			  u3 = u2*u;
			  dwk_j = normkern4 * (4.0*Afac*u3 + 3.0*Bfac*u2 + 2.0*Cfac*u + Dfac);
			}
		    }
                  else
                    {
                      dwk_j = 0; 
                    }
#endif

#ifdef WENDLAND
                  if(r2 < h_i2)
                    {
                      hinv = 1.0 / h_i;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
#ifdef WC2
		      fivetim = (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      dwk_i = -hinv4 * normkern4 * fivetim * u * (1.0 + 5.0*u);
#endif
#ifdef WC4
		      fivetim = (1.0 - u) * (1.0 - u) * (1.0 - u);
		      dwk_i = -hinv4 * normkern4 * fivetim * u;
#endif
#ifdef WC6
		      fivetim = (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) *
			(1.0 - u) * (1.0 - u);
		      u2 = u*u;
		      u3 = u2*u;
		      dwk_i = -hinv4 * normkern4 * fivetim *
			(22.0*u + 154.0*u2 + 352.0*u3);
#endif
                    }
                  else
                    {
                      dwk_i = 0;
                    }

                  if(r2 < h_j * h_j)
                    {
                      hinv = 1.0 / h_j;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
#ifdef WC2
		      fivetim = (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      dwk_j = -hinv4 * normkern4 * fivetim * u * (1.0 + 5.0*u);
#endif
#ifdef WC4
		      fivetim = (1.0 - u) * (1.0 - u) * (1.0 - u);
		      dwk_j = -hinv4 * normkern4 * fivetim * u;
#endif
#ifdef WC6
		      fivetim = (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 - u) *
			(1.0 - u) * (1.0 - u);
		      u2 = u*u;
		      u3 = u2*u;
		      dwk_j = -hinv4 * normkern4 * fivetim *
			(22.0*u + 154.0*u2 + 352.0*u3);  
#endif
		    }
                  else
                    {
                      dwk_j = 0;
                    }
#endif

#ifdef HIFAST
                  if(r2 < h_i2)
                    {
                      hinv = 1.0 / h_i;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
                      normkern4 = normfac * hinv4;
                      if(u < bet)
                        {
                          stemp = 1.0 - u;
                          atemp = alp - u;
                          btemp = bet - u;
                          stemp3 = stemp*stemp*stemp;
                          atemp3 = Afac*atemp*atemp*atemp;
                          btemp3 = Bfac*btemp*btemp*btemp;
                          dwk_i = -normkern4 * 4.0 * (stemp3 + atemp3 + btemp3);
                        }
		      else if(u < alp)
                        {
                          stemp = 1.0 - u;
                          atemp = alp - u;
                          stemp3 = stemp*stemp*stemp;
                          atemp3 = Afac*atemp*atemp*atemp;
                          dwk_i = -normkern4 * 4.0 * (stemp3 + atemp3);
                        }
                      else
                        {
                          stemp = 1.0 - u;
                          stemp3 = stemp*stemp*stemp;
                          dwk_i = -normkern4 * 4.0 * stemp3;
                        }
                    }
                  else
                    {
                      dwk_i = 0;
                    }

		  if(r2 < h_j * h_j)
                    {
                      hinv = 1.0 / h_j;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
                      normkern4 = normfac * hinv4;

                      if(u < bet)
                        {
                          stemp = 1.0 - u;
                          atemp = alp - u;
                          btemp = bet - u;
                          stemp3 = stemp*stemp*stemp;
                          atemp3 = Afac*atemp*atemp*atemp;
                          btemp3 = Bfac*btemp*btemp*btemp;
			  dwk_j = -normkern4 * 4.0 * (stemp3 + atemp3 + btemp3);
                        }
                      else if(u < alp)
                        {
                          stemp = 1.0 - u;
                          atemp = alp - u;
                          stemp3 = stemp*stemp*stemp;
                          atemp3 = Afac*atemp*atemp*atemp;
                          dwk_j = -normkern4 * 4.0 * (stemp3 + atemp3);
                        }
		      else
                        {
                          stemp = 1.0 - u;
                          stemp3 = stemp*stemp*stemp;
                          dwk_j = -normkern4 * 4.0 * stemp3;
                        }
                    }
                  else
                    {
                      dwk_j = 0;
                    }
#endif

#ifdef TRUNCGAUSS
		  if(r2 < h_i2)
                    {
                      hinv = 1.0 / h_i;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
		      normkern4 = All.NormKernD * hinv4;

		      temp = exp(-(u*u/(All.Lf*All.Lf)));
		      /* dwk_i = normkern4 * (-2.0 * u / All.Lf / All.Lf * temp); */
#ifdef GFORMONE
		      dwk_i = normkern4 * (-2.0 * u / All.Lf / All.Lf * temp + 2.0 / All.Lf / All.Lf * All.StopKern);
#else
                      dwk_i = normkern4 * (-2.0 * u / All.Lf / All.Lf * (temp - All.StopKern));
#endif
		    }
		  else
		    {
		      dwk_i = 0;
		    }

		  if(r2 < h_j * h_j)
                    {
                      hinv = 1.0 / h_j;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
                      normkern4 = All.NormKernD * hinv4;

                      temp = exp(-(u*u/(All.Lf*All.Lf)));
                      /* dwk_j = normkern4 * (-2.0 * u / All.Lf / All.Lf * temp) */;
#ifdef GFORMONE
		      dwk_j = normkern4 * (-2.0 * u / All.Lf / All.Lf * temp + 2.0 / All.Lf / All.Lf * All.StopKern);
#else
                      dwk_j = normkern4 * (-2.0 * u / All.Lf / All.Lf * (temp - All.StopKern));
#endif
                    }
                  else
                    {
                      dwk_j = 0;
                    }
#endif

#ifdef HIGHSPLINE
		  if(r2 < h_i2)
                    {
                      hinv = 1.0 / h_i;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
		      normkern4 = All.Normspl * hinv4;

		      if(u < All.betspl)
			{
			  dwk_i = normkern4 * (-All.nsspl*pow((All.rspl - u),(All.nsspl-1)) - All.Aspl*All.nsspl*pow((All.alpspl-u),(All.nsspl-1))-
					       All.Bspl*All.nsspl*pow((All.betspl-u),(All.nsspl-1)));
			}
		      else
			{
			  if(u < All.alpspl)
			    {
			      dwk_i = normkern4 * (-All.nsspl*pow((All.rspl - u),(All.nsspl-1)) - All.Aspl*All.nsspl*pow((All.alpspl-u),(All.nsspl-1)));
			    }
			  else
			    {
			      if(u < All.rspl)
				{
				  dwk_i = normkern4 * (-All.nsspl*pow((All.rspl - u),(All.nsspl-1)));
				    }
			      else
				{
				  dwk_i = 0;
				}
			    }
			}
		    }
		  else
		    {
		      dwk_i = 0;
		    }

		  if(r2 < h_j * h_j)
                    {
                      hinv = 1.0 / h_j;
                      hinv4 = hinv * hinv * hinv * hinv;

                      u = r * hinv;
                      normkern4 = All.Normspl * hinv4;

                      if(u < All.betspl)
                        {
                          dwk_j = normkern4 * (-All.nsspl*pow((All.rspl - u),(All.nsspl-1)) - All.Aspl*All.nsspl*pow((All.alpspl-u),(All.nsspl-1))-
					       All.Bspl*All.nsspl*pow((All.betspl-u),(All.nsspl-1)));
                        }
                      else
                        {
                          if(u < All.alpspl)
                            {
                              dwk_j = normkern4 * (-All.nsspl*pow((All.rspl - u),(All.nsspl-1)) - All.Aspl*All.nsspl*pow((All.alpspl-u),(All.nsspl-1)));
                            }
                          else
                            {
                              if(u < All.rspl)
                                {
                                  dwk_j = normkern4 * (-All.nsspl*pow((All.rspl - u),(All.nsspl-1)));
				}
                              else
                                {
                                  dwk_j = 0;
                                }
                            }
			}
		    }
                  else
                    {
		      dwk_j = 0;
                    }
#endif

		  /* Stuff for artificial visc. here */
		  if(soundspeed_i + soundspeed_j > maxSignalVel)
		    maxSignalVel = soundspeed_i + soundspeed_j;

		  /* Set up mean rho for symm. dissipation switches */
                  rho_ij = 0.5 * (rho + rhoj);

#ifndef SPHS
		  if(vdotr2 < 0)	                /* ... artificial viscosity */
		    {
		      mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */

                      vsig = soundspeed_i + soundspeed_j - 3 * mu_ij;

		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;

		      f2 =
			fabs(SphP[j].DivVel) / (fabs(SphP[j].DivVel) + SphP[j].CurlVel +
						0.0001 * soundspeed_j / fac_mu / SphP[j].Hsml);

		      visc = 0.25 * All.ArtBulkViscConst * vsig * (-mu_ij) / rho_ij * (f1 + f2);
                    }
                  else
                    {
                      visc = 0;
                    }
#endif

#ifdef SPHS
		  if(vdotr2 < 0)                        /* ... artificial viscosity */
                    {
		      mu_ij = fac_mu * vdotr2 / r;      /* note: this is negative! */		  
		      vsig = soundspeed_i + soundspeed_j - 3 * mu_ij;

		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;
		  
		      alpv_ij = 0.5 * (alpv + SphP[j].Alpv);
		      visc = 0.5 * alpv_ij * vsig * (-mu_ij) / rho_ij;
		    }
		  else
		    {
		      visc = 0;
		    }
#endif
		  
#ifndef NOVISCOSITYLIMITER
		  /* make sure that viscous acceleration is not too large */
		  dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * All.Timebase_interval;
		  if(dt > 0 && (dwk_i + dwk_j) < 0)
		    {
		      visc = dmin(visc, 0.5 * fac_vsic_fix * vdotr2 /
				  (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
		    }
#endif
		  /* .... end artificial viscosity evaluation */
		  hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;


#ifdef SPHS

#ifndef FULLCONSERVATION
                  /* SPHS momentum equation (c.f. OSPH) */
                  hfc = hfc_visc + P[j].Mass * 0.5 * (dwk_i + dwk_j) / r *
                    (pressurej + pressure)/(rho*rhoj);
#else
                  /* Use GADGET2 original momentum equation, but with other SPHS effects still on */
                  p_over_rho2_j *= SphP[j].DhsmlDensityFactor;
                  hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
#endif

#else
                  /* GADGET2 original momentum equation */
		  p_over_rho2_j *= SphP[j].DhsmlDensityFactor;
                  hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
#endif

#ifdef NOEZERO
		  /* Higher order E0 subtracted equation */
		  if(alpv_ij < 0.4)
		    alpv_ij = 0;
		  hfc = hfc_visc + P[j].Mass * 0.5 * (dwk_i + dwk_j) / r *
		    ((pressurej - pressure)/(rho*rhoj) + alpv_ij*2.0*pressure/(rho*rhoj));
#endif
#ifdef NEGPRES
		  /* This currently ONLY works properly for the CS kernel!! */

		  /* Add negative pressure to regularize ... The magic number 0.447 is
		     W_CS(0.5 h) c.f. Monaghan 2000. The other magic number 0.01 describes
		     the strength of the negative pressure. Note that this pressure
		     increases as the particle separation decreases. It should act only 
		     on the kernel scale to ensure particle regularity. */

		  /* TEST */
		  /* 
		     printf("Force before addition: %g | %g %g %g | %g %g\n",
		     hfc,wk_i * h_i * h_i * h_i,wk_j * h_j * h_j * h_j,
		     pow(0.5 * (wk_i*h_i*h_i*h_i + wk_j*h_j*h_j*h_j) / 0.44705479,4),u,r);
		  */

		  hfc += 0.01 * P[j].Mass * 0.5 * (dwk_i + dwk_j) / r *
                    (pressurej + pressure)/(rho*rhoj) * 
		    pow(0.5 * (wk_i*h_i*h_i*h_i + wk_j*h_j*h_j*h_j) / 0.44705479,4);

		  /* 
		     printf("Force after addition: %g\n",
		     hfc);
		     fflush(stdout);
		  */
#endif


		  /* Acceleration update here */ 
		  acc[0] -= hfc * dx;
		  acc[1] -= hfc * dy;
		  acc[2] -= hfc * dz;

		  /* Entropy update here */
		  dtEntropy += 0.5 * hfc_visc * vdotr2;

#ifdef LIMIT_TIMESTEPS
                  /* If we have found an inactive neighbour of this active particle that
                     has a longer timestep than this active particle, then set the
                     timestep for the inactive particle equal to that of this active
                     particle */
                  if(P[j].Ti_endstep != All.Ti_Current)
                    if(P[j].TimeStep > All.TimeStepFac * timestep)
                      P[j].TimeStep = All.TimeStepFac * timestep;
#endif

#ifdef SPHS
		  /* Set up modified signal velocity & pressure limiter */
#ifdef CONVERGDIFF
		  if(vdotr2 < 0)
		    {
		      mu_ij = fac_mu * vdotr2 / r;      /* note: this is negative! */
		      mvsig = soundspeed_i + soundspeed_j - 3 * mu_ij;
		      pressurelimiter = fabs(pressure - SphP[j].Pressure) / (pressure + SphP[j].Pressure);
		      
		      /* Entropy dissipation */
		      alpe_ij = 0.5 * (alpv + SphP[j].Alpv) * All.ArtBulkDissConst / All.ArtBulkViscConst;

		      /* This is the energy conserving form */
		      entdiss += P[j].Mass / rho_ij * alpe_ij * mvsig *
			(entropy - SphP[j].Entropy*pow(rhoj/rho,GAMMA_MINUS1)) *
			(dwk_i + dwk_j) / 2.0 * pressurelimiter;
		      
		      /* Mass dissipation */
		      alpm_ij = 0.5 * (alpv + SphP[j].Alpv) * All.ArtBulkDissConst / All.ArtBulkViscConst;

		      /* Conservative mass dissipation (with P-limiter) */
		      massdiss += 0.5 * (mass + P[j].Mass) / rho_ij * alpm_ij * mvsig *
			(mass - P[j].Mass) * (dwk_i + dwk_j) / 2.0 * pressurelimiter;
		    }
#endif
#ifdef WADSLEYDIFF
		  /* Wadsley diffusion */
		  mvsig =  (h_i + h_j)/r * sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
		  pressurelimiter = fabs(pressure - SphP[j].Pressure) / (pressure + SphP[j].Pressure);

		  /* Entropy dissipation */
		  alpe_ij = 0.5 * (alpv + SphP[j].Alpv) * All.ArtBulkDissConst / All.ArtBulkViscConst;

		  /* This is the energy conserving form */
		  entdiss += P[j].Mass / rho_ij * alpe_ij * mvsig *
		    (entropy - SphP[j].Entropy*pow(rhoj/rho,GAMMA_MINUS1)) *
		    (dwk_i + dwk_j) / 2.0 * pressurelimiter;

		  /* Mass dissipation */
		  alpm_ij = 0.5 * (alpv + SphP[j].Alpv) * All.ArtBulkDissConst / All.ArtBulkViscConst;

		  /* Conservative mass dissipation (with P-limiter) */
		  massdiss += 0.5 * (mass + P[j].Mass) / rho_ij * alpm_ij * mvsig *
		    (mass - P[j].Mass) * (dwk_i + dwk_j) / 2.0 * pressurelimiter;
#endif	      
#ifdef MODNEGDIFF
		  /* Allow some diffusion for diverging particles */
		  mu_ij = fac_mu * vdotr2 / r;      /* note: this is negative! */
		  mvsig = soundspeed_i + soundspeed_j - 3 * mu_ij;
		  if(mvsig < 0)
		    mvsig = 0;
                  pressurelimiter = fabs(pressure - SphP[j].Pressure) / (pressure + SphP[j].Pressure);

                  /* Entropy dissipation */
                  alpe_ij = 0.5 * (alpv + SphP[j].Alpv) * All.ArtBulkDissConst / All.ArtBulkViscConst;

                  /* This is the energy conserving form */
                  entdiss += P[j].Mass / rho_ij * alpe_ij * mvsig *
                    (entropy - SphP[j].Entropy*pow(rhoj/rho,GAMMA_MINUS1)) *
                    (dwk_i + dwk_j) / 2.0 * pressurelimiter;

                  /* Mass dissipation */
                  alpm_ij = 0.5 * (alpv + SphP[j].Alpv) * All.ArtBulkDissConst / All.ArtBulkViscConst;
#ifdef NOMASSDISS
		  /* For testing only ... */
		  alpm_ij = 0;
#endif
                  /* Conservative mass dissipation (with P-limiter) */
                  massdiss += 0.5 * (mass + P[j].Mass) / rho_ij * alpm_ij * mvsig *
                    (mass - P[j].Mass) * (dwk_i + dwk_j) / 2.0 * pressurelimiter;
#endif
#endif
		  
#ifdef ERRORS
		  /* Error analysis here */ 
#ifdef SPHS
		  phi = rho;
		  phij = rhoj;
#else
		  phi = 1.0;
		  phij = 1.0;
#endif
#ifdef FULLCONSERVATION
		  phi = 1.0;
                  phij = 1.0;
#endif
		  error0x += P[j].Mass / rhoj * (dwk_i * h_i + dwk_j * h_j) / 2.0 / r * dx *
		    ((phi/phij)*(rhoj/rho) + (phij/phi)*(rho/rhoj));
                  error0y += P[j].Mass / rhoj * (dwk_i * h_i + dwk_j * h_j) / 2.0 / r * dy * 
		    ((phi/phij)*(rhoj/rho) + (phij/phi)*(rho/rhoj));
                  error0z += P[j].Mass / rhoj * (dwk_i * h_i + dwk_j * h_j) / 2.0 / r * dz * 
		    ((phi/phij)*(rhoj/rho) + (phij/phi)*(rho/rhoj));
#endif
		}
	    }
	}
    }
  while(startnode >= 0);

  /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].HydroAccel[k] = acc[k];
      SphP[target].DtEntropy = dtEntropy;
      SphP[target].MaxSignalVel = maxSignalVel;

#ifdef ERRORS
      SphP[target].Error0[0] = error0x;
      SphP[target].Error0[1] = error0y; 
      SphP[target].Error0[2] = error0z;
#endif

#ifdef SPHS
      SphP[target].EntDiss = entdiss;
      SphP[target].MassDiss = massdiss;
#endif
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];
      HydroDataResult[target].DtEntropy = dtEntropy;
      HydroDataResult[target].MaxSignalVel = maxSignalVel;

#ifdef ERRORS
      HydroDataResult[target].Error0[0] = error0x;
      HydroDataResult[target].Error0[1] = error0y; 
      HydroDataResult[target].Error0[2] = error0z;
#endif

#ifdef SPHS
      HydroDataResult[target].EntDiss = entdiss;
      HydroDataResult[target].MassDiss = massdiss;
#endif
    }

}


/*! This is a comparison kernel for a sort routine, which is used to group
 *  particles that are going to be exported to the same CPU.
 */
int hydro_compare_key(const void *a, const void *b)
{
  if(((struct hydrodata_in *) a)->Task < (((struct hydrodata_in *) b)->Task))
    return -1;
  if(((struct hydrodata_in *) a)->Task > (((struct hydrodata_in *) b)->Task))
    return +1;
  return 0;
}
