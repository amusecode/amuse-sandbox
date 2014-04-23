#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

/* This routine calculates the accelerations for all particles in a 
   fixed potential based on their positions. The FORCES in 
   CARTESIAN coordinates for such a fixed potential must be known. As an 
   option the potentials may also be computed so that the total system 
   energy can be correctly calculated.
 */

#ifdef FIXEDPOT

void fixed_potential() 
{
  int    i;
  double force;
  double M,a,b,dummy1,dummy2,dummy3;
  double v0,rc,q0,rt;
  double unit_t, unit_v, unit_l, unit_m;

  /* Here is a list of the potential potentials! (forces really). 
     Other values for FixPotential (i.e. fixedpotential==0) will be 
     ignored. You can add your own here too!
  */

  /* For testing purposes only, a constant force in the x-direction */
  if(All.FixedPotential==1)
    {
      /* Parameters */
      force=0.00001;

      /* Force */
      for(i=1;i<=NumPart;i++)   
	{
	  if(P[i].Mass==0)
	    {
	      P[i].FixAccel[0] = force/header.mass[P[i].Type];
	    }
	  else
	    P[i].FixAccel[0] = force/P[i].Mass;
	}
    }

  /* A Miyamoto-Nagai disc centred on the origin (see B&T p. 43) */
  if(All.FixedPotential==2)
    {
      /* Parameters */
      unit_t = 0.977453;              /* Gyrs */
      unit_v = 1;                     /* kpc */
      unit_l = 1;                     /* km/s */
      unit_m = 2.33e5;                /* Solar masses */
   
      M=5e10/unit_m;                  /* Mass of disc */
      a=4/unit_l;                     /* Scale length of disc */
      b=0.5/unit_l;                   /* Scale height of disc */

      for(i=1;i<=NumPart;i++)   
	{
	  /* Force */
	  dummy1 = a+sqrt(P[i].Pos[2]*P[i].Pos[2]+b*b);
	  dummy2 = sqrt(P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]
			+dummy1*dummy1);
	  dummy3 = pow(dummy2,3);
	  
	  P[i].FixAccel[0] = -All.G*M*P[i].Pos[0]/dummy3;
	  P[i].FixAccel[1] = -All.G*M*P[i].Pos[1]/dummy3;
	  P[i].FixAccel[2] = -All.G*M*P[i].Pos[2]*dummy1/
	    (dummy3*sqrt(P[i].Pos[2]*P[i].Pos[2]+b*b));

	  /* Potential */
	  P[i].FixPot = -All.G*M/dummy2;
	}
    }

  /* A Miyamoto-Nagai disc PLUS log-halo centred on the origin 
     (see B&T ps. 43 and 46). This is chosen to mimic Sverre's code. */
  if(All.FixedPotential==3)
    {
      /* Parameters */
      unit_t = 0.977453;              /* Gyrs */
      unit_v = 1;                     /* kpc */
      unit_l = 1;                     /* km/s */
      unit_m = 2.33e5;                /* Solar masses */

      M=5e10/unit_m;                  /* Mass of disc */
      a=4/unit_l;                     /* Scale length of disc */
      b=0.5/unit_l;                   /* Scale height of disc */

      v0=220/unit_v;                  /* Amplitude of total rotation curve */
      rt=8/unit_l;                    /* defined at infinity and at rt*/
      q0=1;                           /* Halo flattening (1=spherical) */

      /* Calculate log-halo parameters (rc) from these */
      dummy1 = sqrt(rt*rt+(a+b)*(a+b));
      dummy1 = pow(dummy1,3);
      rc = rt*sqrt(v0/sqrt(v0*v0-(rt*rt*All.G*M/dummy1))-1);

      for(i=1;i<=NumPart;i++)   
	{
	  /* Disc Force */
	  dummy1 = a+sqrt(P[i].Pos[2]*P[i].Pos[2]+b*b);
	  dummy2 = sqrt(P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]
			+dummy1*dummy1);
	  dummy3 = pow(dummy2,3);
	  
	  P[i].FixAccel[0] = -All.G*M*P[i].Pos[0]/dummy3;
	  P[i].FixAccel[1] = -All.G*M*P[i].Pos[1]/dummy3;
	  P[i].FixAccel[2] = -All.G*M*P[i].Pos[2]*dummy1/
	    (dummy3*sqrt(P[i].Pos[2]*P[i].Pos[2]+b*b));

	  /* Disc Potential */
	  P[i].FixPot = -All.G*M/dummy2;

	  /* Halo Force */
	  dummy1 = rc*rc+P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+
	    (P[i].Pos[2]*P[i].Pos[2])/(q0*q0);
	  
	  P[i].FixAccel[0] = P[i].FixAccel[0]-v0*v0*P[i].Pos[0]/dummy1;
	  P[i].FixAccel[1] = P[i].FixAccel[1]-v0*v0*P[i].Pos[1]/dummy1;
	  P[i].FixAccel[2] = P[i].FixAccel[2]-v0*v0*P[i].Pos[2]/
	    (dummy1*q0*q0);

	  /* Halo Potential */
	  P[i].FixPot = P[i].FixPot+0.5*v0*v0*log(dummy1);    
	}
    }

  /* A Hernquist alpha profile centred on the origin */
  if(All.FixedPotential==4)
    {
      /* Parameters */
      unit_t = 0.977453;              /* Gyrs */
      unit_v = 1;                     /* kpc */
      unit_l = 1;                     /* km/s */
      unit_m = 2.33e5;                /* Solar masses */
      
      M=1e12/unit_m;                  /* Mass of halo */
      a=25/unit_l;                    /* Scale length of halo */
      b=1;                            /* Central log-slope */

      for(i=1;i<=NumPart;i++)   
	{
	  /* Force */
	  P[i].FixAccel[0] = M*All.G/(2.0-b)*
	    pow(1.0+1/(sqrt(P[i].Pos[0]*P[i].Pos[0]+
			    P[i].Pos[1]*P[i].Pos[1]+
			    P[i].Pos[2]*P[i].Pos[2]))*a,1.0*b-2.0)*
	    (b-2.0)/sqrt(pow(P[i].Pos[0]*P[i].Pos[0]+
			     P[i].Pos[1]*P[i].Pos[1]+
			     P[i].Pos[2]*P[i].Pos[2],3.0))*P[i].Pos[0]/
	    (1.0+1/(sqrt(P[i].Pos[0]*P[i].Pos[0]+
			 P[i].Pos[1]*P[i].Pos[1]+
			 P[i].Pos[2]*P[i].Pos[2]))*a);
	  P[i].FixAccel[1] = M*All.G/(2.0-b)*
	    pow(1.0+1/(sqrt(P[i].Pos[0]*P[i].Pos[0]+
			    P[i].Pos[1]*P[i].Pos[1]+
			    P[i].Pos[2]*P[i].Pos[2]))*a,1.0*b-2.0)*
	    (b-2.0)/sqrt(pow(P[i].Pos[0]*P[i].Pos[0]+
			     P[i].Pos[1]*P[i].Pos[1]+
			     P[i].Pos[2]*P[i].Pos[2],3.0))*P[i].Pos[1]/
	    (1.0+1/(sqrt(P[i].Pos[0]*P[i].Pos[0]+
			 P[i].Pos[1]*P[i].Pos[1]+
			 P[i].Pos[2]*P[i].Pos[2]))*a);
	  P[i].FixAccel[2] = M*All.G/(2.0-b)*
	    pow(1.0+1/(sqrt(P[i].Pos[0]*P[i].Pos[0]+
			    P[i].Pos[1]*P[i].Pos[1]+
			    P[i].Pos[2]*P[i].Pos[2]))*a,1.0*b-2.0)*
	    (b-2.0)/sqrt(pow(P[i].Pos[0]*P[i].Pos[0]+
			     P[i].Pos[1]*P[i].Pos[1]+
			     P[i].Pos[2]*P[i].Pos[2],3.0))*P[i].Pos[2]/
	    (1.0+1/(sqrt(P[i].Pos[0]*P[i].Pos[0]+
			 P[i].Pos[1]*P[i].Pos[1]+
			 P[i].Pos[2]*P[i].Pos[2]))*a);

	  /* Potential */
	  P[i].FixPot = M*All.G/a/(2.0-b)*
	    (pow(1.0+1/(sqrt(P[i].Pos[0]*P[i].Pos[0]+
			     P[i].Pos[1]*P[i].Pos[1]+
			     P[i].Pos[2]*P[i].Pos[2]))*a,1.0*b-2.0)-1.0);
	}
    }
}

#endif
