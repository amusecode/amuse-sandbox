#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file init.c
 *  \brief Code for initialisation of a simulation from initial conditions
 */

#ifdef HIGHSPLINE
/* Short function used by highspline to calculate the kernel norm */
double intsplfunc(double a, double b, double r, int n)
{
  double func;

  func = -(b*b)*pow((r-b),(n+1))/(n+1) - 2.0*b*pow((r-b),(n+2))/(n+1.0)/(n+2.0) - 
    2.0*pow((r-b),(n+3))/(n+1.0)/(n+2.0)/(n+3.0) + 
    (a*a)*pow((r-a),(n+1))/(n+1.0) + 2.0*a*pow((r-a),(n+2))/(n+1.0)/(n+2.0) + 
    2.0*pow((r-a),(n+3))/(n+1.0)/(n+2.0)/(n+3.0);
  
  return(func);
}
#endif

/*! This function reads the initial conditions, and allocates storage for the
 *  tree. Various variables of the particle data are initialised and An intial
 *  domain decomposition is performed. If SPH particles are present, the inial
 *  SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3;

  All.Time = All.TimeBegin;

  switch (All.ICFormat)
    {
    case 1:
#if (MAKEGLASS > 1)
      seed_glass();
#else
      read_ic(All.InitCondFile);
#endif
      break;
    case 2:
    case 3:
      read_ic(All.InitCondFile);
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;
  All.Ti_Current = 0;

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      a3 = All.Time * All.Time * All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      a3 = 1;
    }

  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }

  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;

#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif

#ifdef FIXEDPOT
      for(j = 0; j < 3; j++)
	P[i].FixAccel[j] = 0;
      P[i].FixPot = 0;
#endif

      P[i].Ti_endstep = 0;
      P[i].Ti_begstep = 0;

      P[i].OldAcc = 0;
      P[i].GravCost = 1;
      P[i].Potential = 0;
    }

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

#ifdef FLEXSTEPS
  All.PresentMinStep = TIMEBASE;
  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
    }
#endif


  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].HydroAccel[j] = 0;
#ifdef ERRORS
	  SphP[i].Error0[j] = 0;
#endif
	}

#ifdef SPHS
      SphP[i].EntDiss = 0;
      SphP[i].MassDiss = 0;
#ifndef SPHSAIN
      SphP[i].Alpv = 0;
#endif
      for(j = 0; j < 10; j++)
        {
	  SphP[i].QVx[j] = 0;
          SphP[i].QVy[j] = 0;
          SphP[i].QVz[j] = 0;
	  SphP[i].M0[j] = 0;
        }
      for(j = 0; j < 6; j++)
        SphP[i].M1[j] = 0;
      SphP[i].M2[0] = 0;
      SphP[i].M2[1] = 0;
      SphP[i].M2[2] = 0;
      SphP[i].M3 = 0;
      for(j = 0; j < 6; j++)
        SphP[i].M4[j] = 0;
      for(j = 0; j < 5; j++)
        SphP[i].M5[j] = 0;
      for(j = 0; j < 4; j++)
        SphP[i].M6[j] = 0;
      SphP[i].MaxSignalVel = 0;
      SphP[i].DtAlpv = 0; 
      SphP[i].Alpvloc = 0;

#ifdef SHRINKSTEP
      SphP[i].Hmean = 0;
      SphP[i].Norm = 0;
#endif

#ifdef RHOTHERM
      SphP[i].DenTherm = 0;
#endif

#ifdef SHRINKCOURANT
      SphP[i].CourantMin = 0;
#endif

#endif

#ifdef MASSIFY
      /* For smooth IC generation (using mass noise) */
      SphP[i].DenStore = SphP[i].Density;
#endif

      SphP[i].DtEntropy = 0;

      if(RestartFlag == 0)
	{
	  SphP[i].Hsml = 0;
	  SphP[i].Density = -1;
	}
    }

#ifdef SPHS
  if(ThisTask == 0)
    {
      printf("\nAllocating memory for fluid interp matrices (init loop) ... \n");
    }
  /* Allocate memory for fluid interp matrices */
#ifndef TWODIMS
  All.vecin = gsl_vector_alloc (10);
  All.vecout = gsl_vector_alloc (10);
  All.matrixlu = gsl_matrix_alloc (10,10);
  All.perm = gsl_permutation_alloc (10);
#else
  All.vecin = gsl_vector_alloc (6);
  All.vecout = gsl_vector_alloc (6);
  All.matrixlu = gsl_matrix_alloc (6,6);
  All.perm = gsl_permutation_alloc (6);
#endif
#endif

#ifdef TRUNCGAUSS
  /* Initialize kernal normalisation and factors [Nmin=32; top of pow] */
  All.Nmin = 32.0;
  All.Lf = pow(All.Nmin/All.DesNumNgb,1.0/3.0);
  All.StopKern = exp(-1.0/(All.Lf*All.Lf));
  All.NormKern = 1.0/(-(All.Lf*All.Lf)*2.0*M_PI*All.StopKern + 
		      All.Lf*All.Lf*All.Lf*pow(M_PI,3.0/2.0)*erf(1.0/All.Lf)-4.0/3.0*M_PI*All.StopKern);

#ifdef GFORMONE
  All.NormKernD = -1.0/(pow(M_PI,3.0/2.0)*(All.Lf*All.StopKern + All.Lf*All.Lf*All.Lf*(All.StopKern-1.0)+
					   2.0/3.0/All.Lf*All.StopKern));
#else
  All.NormKernD = -1.0/(pow(M_PI,3.0/2.0)*(All.Lf*All.StopKern + All.Lf*All.Lf*All.Lf*(All.StopKern-1.0)+
                                           1.0/2.0/All.Lf*All.StopKern));
#endif

  /* TEST : Output coefficients: */
  if(ThisTask == 0)
    {
      printf("\n\nTruncated Gaussian coefficients for neighbour number: %g\n",All.DesNumNgb);
      printf("Nmin: %g, Lf: %g, StopKern: %g, NormKern: %g, NormKernD: %g\n\n",All.Nmin,All.Lf,
	     All.StopKern,All.NormKern,All.NormKernD);
      fflush(stdout);
    }
#endif

#ifdef HIGHSPLINE
  /* Initialise highspline kernel */
  All.rspl = 1.0;
  All.Nmin = 64.0;
  All.betspl = pow(All.Nmin/All.DesNumNgb,1.0/3.0);
  All.alpspl = All.betspl + 3.0 * (All.rspl - All.betspl) / 4.0;
  All.nsspl = 2 + round(pow(All.DesNumNgb / All.Nmin,1.0/3.0));
  All.Aspl = pow(All.rspl,All.nsspl-3) * (All.rspl*All.rspl - All.betspl*All.betspl) / 
    (pow(All.alpspl,All.nsspl-3)*(All.alpspl*All.alpspl - All.betspl*All.betspl));
  All.Bspl = -(pow(All.rspl,All.nsspl-1) + All.Aspl*pow(All.alpspl,All.nsspl-1)) / (pow(All.betspl,All.nsspl-1));
  All.Normspl = 1.0/(4.0*M_PI*(intsplfunc(0.0,All.betspl,All.rspl,All.nsspl) + All.Aspl*intsplfunc(0.0,All.betspl,All.alpspl,All.nsspl)+
			       All.Bspl*intsplfunc(0.0,All.betspl,All.betspl,All.nsspl)+intsplfunc(All.betspl,All.alpspl,All.rspl,All.nsspl)+
			       All.Aspl*intsplfunc(All.betspl,All.alpspl,All.alpspl,All.nsspl) + intsplfunc(All.alpspl,All.rspl,All.rspl,All.nsspl)));

  /* TEST : Output coefficients: */
  if(ThisTask == 0)
    {
      printf("\n\nHigh spline coefficients for neighbour number: %g\n",All.DesNumNgb);
      printf("Nmin: %g, Beta: %g, Alpha: %g, Ns: %d\n",All.Nmin,All.betspl,All.alpspl,All.nsspl);
      printf("A: %g, B: %g, Norm: %g\n\n",All.Aspl,All.Bspl,All.Normspl);
      fflush(stdout);
    }
#endif

  ngb_treeallocate(MAX_NGB);

  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  ngb_treebuild();		/* will build tree */

  setup_smoothinglengths();

  TreeReconstructFlag = 1;

  /* at this point, the entropy variable normally contains the 
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly. 
   * Once the density has been computed, we can convert thermal energy to entropy.
   */
#ifndef ISOTHERM_EQS
#ifndef ENTIN

  if(header.flag_entropy_instead_u == 0)
    {
      for(i = 0; i < N_gas; i++)
	{
	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / 
	    pow(SphP[i].Density / a3, GAMMA_MINUS1);
	}
    }

#endif
#endif

}


/*! This routine computes the mass content of the box and compares it to the
 *  specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

#ifdef SPHS
#ifdef ICSMOOTH
  int niter, cnt;
  double psum;
  int iter;
#endif
#endif
#ifdef MASSIFY
  int niter, cnt;
  double Mmin, Mmax, temp;
#endif
#ifdef GLASSIFY
  int niter, cnt, k;
#endif

  if(RestartFlag == 0)
    {

      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }
#ifndef TWODIMS
	  SphP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  SphP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
	}
    }

  if(ThisTask == 0)
    {
      printf("Doing very first density smooth ... \n");
      fflush(stdout);
    }
  density();

#ifdef MASSIFY
  /* Iterate the particle masses to give const. (or desired) density */
  niter = 50;
  Mmin = 1e30;
  Mmax = 0;
  for(i = 0; i < N_gas; i++)
    {
      if(P[i].Mass > Mmax)
	Mmax = P[i].Mass;
      if(P[i].Mass < Mmin)
	Mmin = P[i].Mass;
    }

  if(ThisTask == 0)
    {
      printf("*** WARNING *** SPHS MASSIFY density iterate ... \n");
      printf(" ... I hope this is what you want ... \n");
      fflush(stdout);
    }
  
  for(cnt = 0; cnt < niter; cnt++)
    {
      for(i = 0; i < N_gas; i++)
	{
	  temp = P[i].Mass * SphP[i].DenStore / SphP[i].Density;
	  if(temp < 5.0 * Mmax && temp > Mmin / 5.0)
	    P[i].Mass = temp;
	}
      density();
      printf("Iteration: %d; Target: %g; True: %g\n",cnt,SphP[i-1].Density,SphP[i-1].DenStore);
      fflush(stdout);
    }
#endif

#ifdef GLASSIFY
  /* WARNING: THIS DOES NOT WORK YET IN THIS VERSION !!! */


  /* Generate glass IC from analytic density profile (in densityfunc; hydra.c) */
  niter = 1500;
  
  if(ThisTask == 0)
    {
      printf("*** WARNING *** SPHS GLASSIFY iterate ... \n");
      printf(" ... I hope this is what you want ... \n");
      fflush(stdout);
    }
  
  for(cnt = 0; cnt < niter; cnt++)
    {
      /* Calcualte the density profile and drift forces */
      density();
      force_update_hmax();
      hydra();
      
      /* Drift particles */
      for(i = 0; i < N_gas; i++)
	for(k = 0; k < 3; k++)
	  P[i].Pos[k] += All.GlassTol * SphP[i].HydroAccel[k];

      /* Some output */
      if(ThisTask == 0)
	{
	  printf("Iteration: %d | Density: %g | Accel: (%g,%g,%g)\n",cnt,SphP[i-1].Density, SphP[i-1].HydroAccel[0],SphP[i-1].HydroAccel[1],SphP[i-1].HydroAccel[2]);
	  fflush(stdout);
	}
    }
#endif

#ifdef ICSMOOTH
#ifdef SPHS
  /* Set the entropy to give pressure equilibrium exactly initially */
  /* Initial constant pressure in sim units. This is set in the parameter file in 
     FixedPressure (double) */
  if(All.FixedPressure > 0)
    {
      if(ThisTask == 0)
        {
          printf("*** WARNING *** IC smoothing with constant pressure: %g.\n",
                 All.FixedPressure);
          printf(" ... I hope this is what you want ... \n");
          fflush(stdout);
        }

      for(i = 0; i < N_gas; i++)
        {
          SphP[i].Entropy = All.FixedPressure / pow(SphP[i].Density,GAMMA);
          SphP[i].Pressure = All.FixedPressure;
	}

#ifdef RHOTHERM
      /* Iterate to obtain constant pressure [if RHOTHERM on] */
      for (iter = 0; iter < 5; iter++)
	{
	  if(ThisTask == 0)
	    {
	      printf("Pressure smoothing iteration %d | %g ... \n",iter,psum);
	      fflush(stdout);
	    }
	  
	  density();
	  
	  psum = 0;
	  for(i = 0; i < N_gas; i++)
	    {
	      psum += SphP[i].Entropy * pow(SphP[i].DenTherm,GAMMA);
	      SphP[i].Entropy = All.FixedPressure / pow(SphP[i].DenTherm,GAMMA);
	    }
	  psum = psum / (N_gas * All.FixedPressure);
	}
#endif
    }
#endif
#endif

#ifdef SPHS
  force_update_hmax(); 

  /* Need to run twice the very first time to calculate
     the signal velocity. From then on, we use the
     value from the previous timestep */
  hydro_force();
  density();
  hydro_force();

  /* 
     For this first interpolation (and only this one) set 
     all dissipation terms to maximum if there is a discontinuity
     in any fluid variable. This is to cope with situations
     like the Sod shock tube that has a pressure discontinuity 
     but no initial velocity discontinuity. If the visc is not
     required at the discontinuity then it will rapidly decay
     away. Depending on your problem, you may want to change
     this. 
  */

#ifdef ICALPHASOD
  for(i = 0; i < N_gas; i++)
    {
      if(P[i].Pos[2] > 0.95 && P[i].Pos[2] < 1.05)
	{
	  SphP[i].Alpv = All.ArtBulkViscConst;
	  SphP[i].Alpvloc = All.ArtBulkViscConst;
	}
    }
#endif
#endif

}


/*! If the code is run in glass-making mode, this function populates the
 *  simulation box with a Poisson sample of particles.
 */
#if (MAKEGLASS > 1)
void seed_glass(void)
{
  int i, k, n_for_this_task;
  double Range[3], LowerBound[3];
  double drandom, partmass;
  long long IDstart;

  All.TotNumPart = MAKEGLASS;
  partmass = All.Omega0 * (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G))
    * (All.BoxSize * All.BoxSize * All.BoxSize) / All.TotNumPart;

  All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);	/* sets the maximum number of particles that may */

  allocate_memory();

  header.npartTotal[1] = All.TotNumPart;
  header.mass[1] = partmass;

  if(ThisTask == 0)
    {
      printf("\nGlass initialising\nPartMass= %g\n", partmass);
      printf("TotNumPart= %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
    }

  /* set the number of particles assigned locally to this task */
  n_for_this_task = All.TotNumPart / NTask;

  if(ThisTask == NTask - 1)
    n_for_this_task = All.TotNumPart - (NTask - 1) * n_for_this_task;

  NumPart = 0;
  IDstart = 1 + (All.TotNumPart / NTask) * ThisTask;

  /* split the temporal domain into Ntask slabs in z-direction */

  Range[0] = Range[1] = All.BoxSize;
  Range[2] = All.BoxSize / NTask;
  LowerBound[0] = LowerBound[1] = 0;
  LowerBound[2] = ThisTask * Range[2];

  srand48(ThisTask);

  for(i = 0; i < n_for_this_task; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  drandom = drand48();

	  P[i].Pos[k] = LowerBound[k] + Range[k] * drandom;
	  P[i].Vel[k] = 0;
	}

      P[i].Mass = partmass;
      P[i].Type = 1;
      P[i].ID = IDstart + i;

      NumPart++;
    }
}
#endif
