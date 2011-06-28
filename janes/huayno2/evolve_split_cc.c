/*
 * Connected components Hamiltonian split (CC-split)
 */

//#define CC2_SPLIT_CONSISTENCY_CHECKS
//#define CC2_SPLIT_SHORTCUTS

#define IS_ZEROSYS(SYS) (((SYS)->n == 0) && ((SYS)->part == NULL) && ((SYS)->last == NULL) && ((SYS)->next_cc == NULL))
#define IS_ZEROSYSs(SYS) (((SYS).n == 0) && ((SYS).part == NULL) && ((SYS).last == NULL) && ((SYS).next_cc == NULL))
#define LOG_CC_SPLIT(C, R) { \
	LOG("clevel = %d s.n = %d c.n = {", clevel, s.n); \
	for (struct sys *_ci = (C); !IS_ZEROSYS(_ci); _ci = _ci->next_cc) printf(" %d ", _ci->n ); \
	printf("} r.n = %d\n", (R)->n); \
};

#define LOGSYS_ID(SYS) for (UINT i = 0; i < (SYS).n; i++) { printf("%u ", (SYS).part[i].id); } printf("\n");
#define LOGSYSp_ID(SYS) for (UINT i = 0; i < (SYS)->n; i++) { printf("%u ", (SYS)->part[i].id); } printf("\n");
#define LOGSYSC_ID(SYS) for (struct sys *_ci = &(SYS); !IS_ZEROSYS(_ci); _ci = _ci->next_cc) {printf("{"); for (UINT i = 0; i < _ci->n; i++) {printf("%u ", _ci->part[i].id); } printf("}\t");} printf("\n");

DOUBLE timestep_ij(struct sys r, UINT i, struct sys s, UINT j) {

	FLOAT timestep;
	FLOAT dx[3],dr3,dr2,dr,dv[3],dv2,mu,vdotdr2,tau,dtau;

	//if(s.part[i].timestep !=HUGE_VAL) endrun((char *)"timestep??");

    timestep=HUGE_VAL;

    dx[0]=r.part[i].pos[0]-s.part[j].pos[0];
    dx[1]=r.part[i].pos[1]-s.part[j].pos[1];
    dx[2]=r.part[i].pos[2]-s.part[j].pos[2];
    dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;

    if(dr2>0) {
    	dr=sqrt(dr2);
    	dr3=dr*dr2;
        dv[0]=r.part[i].vel[0]-s.part[j].vel[0];
        dv[1]=r.part[i].vel[1]-s.part[j].vel[1];
        dv[2]=r.part[i].vel[2]-s.part[j].vel[2];
        vdotdr2=(dv[0]*dx[0]+dv[1]*dx[1]+dv[2]*dx[2])/dr2;
        dv2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
        mu=r.part[i].mass+s.part[j].mass;

#ifdef RATIMESTEP
        tau=RARVRATIO*dt_param/M_SQRT2*sqrt(dr3/mu);
        dtau=3/2.*tau*vdotdr2;
        if(dtau>1.) dtau=1.;
        tau/=(1-dtau/2);
        if(tau < timestep) timestep=tau;
#endif

#ifdef RVTIMESTEP
        if(dv2>0) {
        	tau=dt_param*dr/sqrt(dv2);
        	dtau=tau*vdotdr2*(1+mu/(dv2*dr));
        	if(dtau>1.) dtau=1.;
        	tau/=(1-dtau/2);
        	if(tau < timestep) timestep=tau;
        }
#endif
    }

    if (timestep < 0) {
		endrun((char*)"timestep_ij: negative timestep!\n");
	}

    return timestep;
}

DOUBLE timestep_ij_bw(struct sys r, UINT i, struct sys s, UINT j) {

	FLOAT timestep;
	FLOAT dx[3],dr3,dr2,dr,dv[3],dv2,mu,vdotdr2,tau,dtau;

	//if(s.part[i].timestep !=HUGE_VAL) endrun((char *)"timestep??");

    timestep=HUGE_VAL;

    dx[0]=r.part[i].pos[0]-s.part[j].pos[0];
    dx[1]=r.part[i].pos[1]-s.part[j].pos[1];
    dx[2]=r.part[i].pos[2]-s.part[j].pos[2];
    dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;

    if(dr2>0) {
    	dr=sqrt(dr2);
    	dr3=dr*dr2;
        dv[0]=r.part[i].vel[0]-s.part[j].vel[0];
        dv[1]=r.part[i].vel[1]-s.part[j].vel[1];
        dv[2]=r.part[i].vel[2]-s.part[j].vel[2];
        vdotdr2=(dv[0]*dx[0]+dv[1]*dx[1]+dv[2]*dx[2])/dr2;
        dv2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
        mu=r.part[i].mass+s.part[j].mass;

#ifdef RATIMESTEP
        tau=RARVRATIO*dt_param/M_SQRT2*sqrt(dr3/mu);
        dtau=3/2.*tau*vdotdr2;
        if(dtau<-1.) dtau=-1.;
        tau/=(1 + dtau/2);
        if(tau < timestep) timestep=tau;
#endif

#ifdef RVTIMESTEP
        if(dv2>0) {
        	tau=dt_param*dr/sqrt(dv2);
        	dtau=tau*vdotdr2*(1+mu/(dv2*dr));
        	if(dtau<-1.) dtau=-1.;
        	tau/=(1+dtau/2);
        	if(tau < timestep) timestep=tau;
        }
#endif
    }

    if (timestep < 0) {
		endrun((char*)"timestep_ij_bw: negative timestep!\n");
	}

    return -timestep;
}


void split_cc(struct sys s, struct sys *c, struct sys *r, DOUBLE dt) {
/*
	find_cc: finds a non-trivial connected component, or sets c to NULLSYS, starting from i
	input: system s, pointers+allocated memory for c and r
	output: c is a connected component, r
	particles in s are re-arranged: c is in the beginning, r is at the end
	c end signal is when c->next_cc == zerosys

	Complexity of the algorithm is n**2, but it does not cache time steps;
	Algorithm has to be re-run at each evolve.
*/
	/*
	LOG("Adjacency matrix:\n");
	for (UINT i = 0; i < s.n; i++) {
		LOG("");
		for (UINT j = 0; j < s.n; j++) {
			printf("%f\t", timestep_ij(s, i, j));
		}
		printf("\n");
	}
	*/
	tstep[clevel]++; // not directly comparable to corresponding SF-split statistics

	struct sys *c_next;
	c_next = c;
	*c_next = zerosys;

	UINT processed = 0; // increase if something is added from the stack to the cc

	UINT comp_next = 0; // increase to move a particle from stack to cc; points to the first element of the stack
	UINT comp_size = 0; // amount of particles added to the current cc
	UINT stack_next = 1; // swap this with s[i] to increase the stack
	UINT stack_size = 1; // first element of the stack is s[comp_next]
	                     // last element of the stack is s[comp_next + stack_size - 1]

	UINT rest_next = s.n - 1; // swap this to add to the rest-system

	// find all connected components
	while (processed < s.n) {
		//LOG("split_cc: searching for connected components: %d / %d\n", processed, s.n);

		// search for the next connected component
		while (stack_size > 0) {

			// iterate over all unvisited elements
			for (UINT i = stack_next; i <= rest_next; i++) {

				// if element is connected to the first element of the stack
				DOUBLE timestep = (dt > 0) ? timestep_ij(s, comp_next, s, i) : timestep_ij_bw(s, comp_next, s, i);
				//LOG("timestep: %LE %LE %LE\n", timestep, timestep_ij(s, comp_next, s, i), timestep_ij_bw(s, comp_next, s, i));
				tcount[clevel]++;

				if (((dt > 0) && (timestep <= dt)) ||
				   ((dt < 0) && (timestep >= dt))) {
					//LOG("swap: %d<->%d\n", comp_next, i);

					// add i to the end of the stack by swapping stack_next and i
					SWAP( s.part[ stack_next ], s.part[i], struct particle );
					stack_next++;
					stack_size++;
				}
			}

			// pop the stack; add to the connected component
			comp_size++;
			comp_next++;
			stack_size--;
		}

		processed += comp_size;

		//LOG("split_cc: found component with size: %d\n", comp_size);
		// create a new component with size
		if (comp_size > 1) {
			//LOG("found non-trivial component with size: %d\n", comp_size);
			// create new component c from u[0] to u[cc_visited - 1]
			// remove components from u (u.n, u.part)

			c_next->n = comp_size;
			c_next->part = &( s.part[ comp_next - comp_size ]);
			c_next->last = &( s.part[ comp_next - 1 ]);
			c_next->next_cc = (struct sys*) malloc( sizeof(struct sys) );
			c_next = c_next->next_cc;
			*c_next = zerosys;

			comp_next = stack_next;
			comp_size = 0;
			stack_next = stack_next + 1;
			stack_size = 1;

		// component is trivial: add to rest
		} else {
			//LOG("found trivial component; adding to rest\n");
			SWAP(s.part[ comp_next - 1 ], s.part[ rest_next ], struct particle );
			rest_next--;

			comp_next = comp_next - 1;
			comp_size = 0;
			stack_next = comp_next + 1;
			stack_size = 1;
		}
	}

	if (processed != s.n) {
		LOG("split_cc: processed=%u s.n=%u r->n=%u\n", processed, s.n, r->n);
		endrun((char*)"split_cc: zomg epic fail!\n");
	}

	// create r system
	r->n = (s.n - 1) - rest_next;
	if (r->n > 0) {
		r->part = &( s.part[rest_next + 1] );
		r->last = s.last;
	} else {
		r->part = NULL;
		r->last = NULL;
	}
	//LOG("split_cc: rest system size: %d\n", r->n);

}


void split_cc_verify(struct sys s, struct sys *c, struct sys *r) {
	//LOG("split_cc_verify ping s.n=%d r->n=%d\n", s.n, r->n);
	//LOG_CC_SPLIT(c, r);

	UINT pcount_check = 0;
	for (UINT i = 0; i < s.n; i++) {

		pcount_check = 0;
		UINT particle_found = 0;

		struct particle *p = &( s.part[i] );

	    for (struct sys *cj = c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
			pcount_check += cj->n;
			//LOG("%d\n", pcount_check);

			// search for p in connected components
			for (UINT k = 0; k < cj->n; k++) {
				struct particle * pk = &( cj->part[k] );

				// is pk equal to p
				if (p->id == pk->id) {
					particle_found += 1;
					//LOG("split_cc_verify: found in a cc\n");
				}
			}

			if (& ( cj->part[cj->n - 1] ) != cj->last) {
				LOG("split_cc_verify: last pointer for c is not set correctly!");
				LOG_CC_SPLIT(c, r);
				endrun((char*)"data structure corrupted\n");
			}
	    }

		// search for p in rest
		for (UINT k = 0; k < r->n; k++) {
			struct particle * pk = &( r->part[k] );

			// is pk equal to p
			if (p->id == pk->id) {
				particle_found += 1;
			}
		}

		if (particle_found != 1) {
			LOG("split_cc_verify: particle %d particle_found=%d ", i, particle_found);
			LOG_CC_SPLIT(c, r);
			endrun((char*)"data structure corrupted\n");
		}

	}

	//if (& ( r->part[r->n - 1] ) != r->last) {
	//	LOG("split_cc_verify: last pointer for r is not set correctly! %d %d",&( r->part[r->n - 1] ), r->last);
	//	LOG_CC_SPLIT(c, r);
	//	endrun((char*)"data structure corrupted\n");
	//}

	if (pcount_check + r->n != s.n) {
		LOG("split_cc_verify: particle count mismatch (%d %d)\n", pcount_check + r->n, s.n);
		LOG_CC_SPLIT(c, r);
		endrun((char*)"data structure corrupted\n");
		//endrun((char*)"split_cc_verify: particle count mismatch\n");
	} else {
		//LOG("split_cc_verify pong\n");
	}

	//endrun((char*)"Fin.\n");

}

void split_cc_verify_ts(struct sys *c, struct sys *r, DOUBLE dt) {

	DOUBLE ts_ij;
	// verify C-C interactions
    for (struct sys *ci = c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
    	for (UINT i = 0; i < ci->n; i++) {
       		for (struct sys *cj = c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
    			if (ci == cj) {
    				continue;
    			}
    	    	for (UINT j = 0; j < cj->n; j++) {
    	    		ts_ij = timestep_ij(*ci, i, *cj, j);
    	    		//LOG("comparing %d %d\n", ci->part[i].id, cj->part[j].id);
    	    		//LOG("%f %f \n", ts_ij, dt);
    	    		if (dt > ts_ij) {
    	    			endrun((char*)"split_cc_verify_ts C-C timestep underflow\n");
    	    		}
        		}
    		}
    	}
    }

	// verify C-R interactions
    for (struct sys *ci = c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
    	for (UINT i = 0; i < ci->n; i++) {

    		for (UINT j = 0; j < r->n; j++) {
	    		ts_ij = timestep_ij(*ci, i, *r, j);
    	   		if (ts_ij < dt) {
    	   			endrun((char*)"split_cc_verify_ts C-R timestep underflow\n");
    	   		}
        	}
    	}
    }

    // verify R interactions
	for (UINT i = 0; i < r->n; i++) {
    	for (UINT j = 0; j < r->n; j++) {
    		if (i == j) continue;
    		ts_ij = timestep_ij(*r, i, *r, j);
    		if (ts_ij < dt) {
    			endrun((char*)"split_cc_verify_ts R-R timestep underflow\n");
    		}
		}
	}

}

void free_sys(struct sys * s) {
	if (s==NULL) return;
	if (s->next_cc != NULL) {
		free_sys(s->next_cc);
	}
	free(s);
}

DOUBLE sys_forces_max_timestep(struct sys s) {
  DOUBLE ts = 0.0;
  DOUBLE ts_ij;
  for (UINT i = 0; i < s.n; i++) {
    for (UINT j = 0; j < s.n; j++) {
      if (i != j) {
        ts_ij = timestep_ij(s, i, s, j);
        if (ts_ij >= ts) { ts = ts_ij; };
      }
    }
  }
  return ts;
}

void evolve_split_cc2(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {

	struct sys c = zerosys, r = zerosys;
	clevel++;
	if (etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small\n");

#ifdef CC2_SPLIT_CONSISTENCY_CHECKS
	if (clevel == 0) {
		printf("consistency_checks: ", s.n, clevel);
	}
#endif

#ifdef CC2_SPLIT_CONSISTENCY_CHECKS
	// debug: make a copy of s to verify that the split has been done properly
	struct sys s_before_split;
	s_before_split.n = s.n;
	s_before_split.part = (struct particle*) malloc(s.n*sizeof(struct particle));
	s_before_split.last = &( s_before_split.part[s.n - 1] );
	s_before_split.next_cc = NULL;
	memcpy(s_before_split.part, s.part, s.n*sizeof(struct particle));
#endif

	split_cc(s, &c, &r, dt);
	//if (s.n != c.n) LOG_CC_SPLIT(&c, &r); // print out non-trivial splits

#ifdef CC2_SPLIT_CONSISTENCY_CHECKS
/*
  	if (s.n != r.n) {
		LOG("s: ");
		LOGSYS_ID(s_before_split);
		LOG("c: ");
		LOGSYSC_ID(c);
		LOG("r: ");
		LOGSYS_ID(r);
	}
*/
	// verify the split
	split_cc_verify(s_before_split, &c, &r);
	split_cc_verify_ts(&c, &r, dt);
	free(s_before_split.part);
	if (clevel == 0) {
		printf("ok ");
	}
#endif

#ifdef CC2_SPLIT_SHORTCUTS
	if (s.n == c.n) {
		DOUBLE initial_timestep = sys_forces_max_timestep(s);
		DOUBLE dt_step = dt;

		while (dt_step > initial_timestep) dt_step = dt_step / 2;

		LOG("CC2_SPLIT_SHORTCUTS clevel=%d dt/dt_step=%Le\n", clevel, dt / dt_step);
		for (DOUBLE dt_now = 0; dt_now < dt; dt_now += dt_step) {
			evolve_split_cc2(s, dt_now, dt_now + dt_step,(DOUBLE) dt_step);
		}

		clevel--;
		free_sys(c.next_cc);
		return;
	}
#endif

	if (IS_ZEROSYSs(c)) {
		deepsteps++;
		simtime+=dt;
	}

	// evolve all fast components, 1st time
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		evolve_split_cc2(*ci, stime, stime+dt/2, dt/2);
	}

	drift(r, stime+dt/2, dt/2); // drift r, 1st time

	// kick ci <-> cj
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		for (struct sys *cj = &c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
			if (ci != cj) {
				kick(*ci, *cj, dt);
				//kick(*cj, *ci, dt);
			}
		}
	}

	// kick c <-> rest
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		kick(r, *ci, dt);
		kick(*ci, r, dt);
	}

	kick(r, r, dt); // kick rest

	drift(r, etime, dt/2); 	// drift r, 2nd time

	// evolve all fast components, 2nd time
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		evolve_split_cc2(*ci, stime+dt/2, etime, dt/2);
	}

	clevel--;
	free_sys(c.next_cc);
}

void evolve_split_cc2_twobody(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {
  // TODO fail/warn if smoothing is non-zero
  // TODO use kepler solver only if dt is larger than the orbital period
  if (s.n == 2) {
    //LOG("evolve: close encounters!\n");
    evolve_twobody(s, stime, etime, dt);
  } else {
    //LOG("evolve: regular!\n");
    struct sys c = zerosys, r = zerosys;
    clevel++;
    if (etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small\n");

#ifdef CC2_SPLIT_CONSISTENCY_CHECKS
    if (clevel == 0) {
      printf("consistency_checks: ", s.n, clevel);
    }
#endif

#ifdef CC2_SPLIT_CONSISTENCY_CHECKS
    // debug: make a copy of s to verify that the split has been done properly
    struct sys s_before_split;
    s_before_split.n = s.n;
    s_before_split.part = (struct particle*) malloc(s.n*sizeof(struct particle));
    s_before_split.last = &( s_before_split.part[s.n - 1] );
    s_before_split.next_cc = NULL;
    memcpy(s_before_split.part, s.part, s.n*sizeof(struct particle));
#endif

    split_cc(s, &c, &r, dt);
    //if (s.n != c.n) LOG_CC_SPLIT(&c, &r); // print out non-trivial splits

#ifdef CC2_SPLIT_CONSISTENCY_CHECKS
  /*
      if (s.n != r.n) {
      LOG("s: ");
      LOGSYS_ID(s_before_split);
      LOG("c: ");
      LOGSYSC_ID(c);
      LOG("r: ");
      LOGSYS_ID(r);
    }
  */
    // verify the split
    split_cc_verify(s_before_split, &c, &r);
    split_cc_verify_ts(&c, &r, dt);
    free(s_before_split.part);
    if (clevel == 0) {
      printf("ok ");
    }
#endif

#ifdef CC2_SPLIT_SHORTCUTS
    if (s.n == c.n) {
      DOUBLE initial_timestep = sys_forces_max_timestep(s);
      DOUBLE dt_step = dt;

      while (dt_step > initial_timestep) dt_step = dt_step / 2;

      LOG("CC2_SPLIT_SHORTCUTS clevel=%d dt/dt_step=%Le\n", clevel, dt / dt_step);
      for (DOUBLE dt_now = 0; dt_now < dt; dt_now += dt_step) {
        evolve_split_cc2(s, dt_now, dt_now + dt_step,(DOUBLE) dt_step);
      }

      clevel--;
      free_sys(c.next_cc);
      return;
    }
#endif

    if (IS_ZEROSYSs(c)) {
      deepsteps++;
      simtime+=dt;
    }

    // evolve all fast components, 1st time
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      evolve_split_cc2_twobody(*ci, stime, stime+dt/2, dt/2);
    }

    drift(r, stime+dt/2, dt/2); // drift r, 1st time

    // kick ci <-> cj
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      for (struct sys *cj = &c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
        if (ci != cj) {
          kick(*ci, *cj, dt);
          //kick(*cj, *ci, dt);
        }
      }
    }

    // kick c <-> rest
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      kick(r, *ci, dt);
      kick(*ci, r, dt);
    }

    kick(r, r, dt); // kick rest

    drift(r, etime, dt/2);  // drift r, 2nd time

    // evolve all fast components, 2nd time
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      evolve_split_cc2_twobody(*ci, stime+dt/2, etime, dt/2);
    }

    clevel--;
    free_sys(c.next_cc);
  }
}

void evolve_split_cc4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
void evolve_split_cc4_chunks(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int chunks) {
  //LOG("\nsys_initial_timestep=%Le, dt=%e\n", initial_timestep, dt);
  DOUBLE dt_step = dt / chunks;
  DOUBLE stime_i = stime;
  for (int i = 0; i < chunks; i++) {
    //LOG("chunk %d\n", chunks);
    evolve_split_cc4(s, stime_i, stime_i + dt_step, (DOUBLE)dt_step);
    stime_i += dt_step;
  }
  //LOG("making: %Le steps\n", dt / dt_step);
}

void evolve_split_cc4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {

	//LOG("evolve_split_hold_dkd_cc: stime=%LE etime=%LE dt=%LE\n", stime, etime, dt);
	struct sys c = zerosys, r = zerosys;
	clevel++;

	if (((dt > 0) && (etime <= stime)) || ((dt < 0) && (etime >= stime)) || // generalization of: etime <= stime
		dt==0 || clevel>=MAXLEVEL) endrun((char *)"evolve_split_cc4: timestep too small");

	split_cc(s, &c, &r, dt);
	if (IS_ZEROSYSs(c)) {
		deepsteps++;
		simtime+=dt;
	}

#define EVOLVE_CC(TIME, DT, CHUNKS) \
  for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) { \
    evolve_split_cc4_chunks(*ci, (TIME), (TIME) + (DT), (DT), (CHUNKS)); \
  } \
  (TIME) += (DT);

#define DRIFT_RI(TIME, DT) \
  drift(r, (TIME) + (DT), (DT)); \
  (TIME) += (DT);

#define KICK_RI(TIME, DT) \
  /* kick ci <-> cj */ \
  for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) { \
    for (struct sys *cj = &c; !IS_ZEROSYS(cj); cj = cj->next_cc) { \
      if (ci != cj) { \
        kick(*ci, *cj, (DT)); \
      } \
    } \
  } \
  /* kick c <-> rest */ \
  for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) { \
    kick(r, *ci, (DT)); \
    kick(*ci, r, (DT)); \
  } \
  /* kick rest */ \
  kick(r, r, (DT)); \
  (TIME) += (DT);

	// time-keeping counters
	DOUBLE time_c = stime;
	DOUBLE time_r = stime;
	DOUBLE time_k = stime;
	DOUBLE time_d = stime;

/*
	// 2-nd order split (useful as a test of the evaluation macros)
	#define EVOLVE_SPLIT_DKD2(EVOLVE_A, EVOLVE_B, TIME_A, TIME_B, DT) \
		EVOLVE_A(TIME_A, 0.5*(DT)); \
		EVOLVE_B(TIME_B, 1.0*(DT)); \
		EVOLVE_A(TIME_A, 0.5*(DT));

	#define EVOLVE_SPLIT_CRI2(EVOLVE_A, EVOLVE_B, TIME_A, TIME_B, DT) \
		EVOLVE_A(TIME_A, 0.5*(DT)); \
		EVOLVE_B(TIME_B, 1.0*(DT)); \
		EVOLVE_A(TIME_A, 0.5*(DT));

	#define EVOLVE_RI(TIME, DT) EVOLVE_SPLIT_DKD2(DRIFT_RI, KICK_RI, time_d, time_k, (DT));
		EVOLVE_SPLIT_CRI2(EVOLVE_CC, EVOLVE_RI, time_c, time_r, dt);
*/

	// 4-th order symmetric operator split
	#define A1 0.0792036964311957
	#define A2 0.353172906049774
	#define A3 -.0420650803577195
	#define A4 (1 - 2*((A1) + (A2) + (A3))) // A4=0.21937695575349958
	#define B1 0.209515106613362 // div4
	#define B2 -.143851773179818 // div3
	#define B3 (0.5 - (B1) - (B2)) // B3=0.43433666656645598 // div9

	#define EVOLVE_SPLIT_DKD4(EVOLVE_A, EVOLVE_B, TIME_A, TIME_B, DT) \
		EVOLVE_A(TIME_A, A1*(DT)); \
		EVOLVE_B(TIME_B, B1*(DT)); \
		EVOLVE_A(TIME_A, A2*(DT)); \
		EVOLVE_B(TIME_B, B2*(DT)); \
		EVOLVE_A(TIME_A, A3*(DT)); \
		EVOLVE_B(TIME_B, B3*(DT)); \
		EVOLVE_A(TIME_A, A4*(DT)); \
		EVOLVE_B(TIME_B, B3*(DT)); \
		EVOLVE_A(TIME_A, A3*(DT)); \
		EVOLVE_B(TIME_B, B2*(DT)); \
		EVOLVE_A(TIME_A, A2*(DT)); \
		EVOLVE_B(TIME_B, B1*(DT)); \
		EVOLVE_A(TIME_A, A1*(DT));

#define CC4_REDUCE_FACTOR 1
#define EVOLVE_SPLIT_CRI4(EVOLVE_A, EVOLVE_B, TIME_A, TIME_B, DT) \
    EVOLVE_A(TIME_A, A1*(DT)); \
    EVOLVE_B(TIME_B, B1*(DT), 4 * CC4_REDUCE_FACTOR); \
    EVOLVE_A(TIME_A, A2*(DT)); \
    EVOLVE_B(TIME_B, B2*(DT), 3 * CC4_REDUCE_FACTOR); \
    EVOLVE_A(TIME_A, A3*(DT)); \
    EVOLVE_B(TIME_B, B3*(DT), 9 * CC4_REDUCE_FACTOR); \
    EVOLVE_A(TIME_A, A4*(DT)); \
    EVOLVE_B(TIME_B, B3*(DT), 9 * CC4_REDUCE_FACTOR); \
    EVOLVE_A(TIME_A, A3*(DT)); \
    EVOLVE_B(TIME_B, B2*(DT), 3 * CC4_REDUCE_FACTOR); \
    EVOLVE_A(TIME_A, A2*(DT)); \
    EVOLVE_B(TIME_B, B1*(DT), 4 * CC4_REDUCE_FACTOR); \
    EVOLVE_A(TIME_A, A1*(DT));
/*
#define EVOLVE_SPLIT_CRI4(EVOLVE_A, EVOLVE_B, TIME_A, TIME_B, DT) \
    EVOLVE_A(TIME_A, A1*(DT)); \
    EVOLVE_B(TIME_B, B1*(DT), 1); \
    EVOLVE_A(TIME_A, A2*(DT)); \
    EVOLVE_B(TIME_B, B2*(DT), 1); \
    EVOLVE_A(TIME_A, A3*(DT)); \
    EVOLVE_B(TIME_B, B3*(DT), 1); \
    EVOLVE_A(TIME_A, A4*(DT)); \
    EVOLVE_B(TIME_B, B3*(DT), 1); \
    EVOLVE_A(TIME_A, A3*(DT)); \
    EVOLVE_B(TIME_B, B2*(DT), 1); \
    EVOLVE_A(TIME_A, A2*(DT)); \
    EVOLVE_B(TIME_B, B1*(DT), 1); \
    EVOLVE_A(TIME_A, A1*(DT));
*/
#define EVOLVE_RI(TIME, DT) EVOLVE_SPLIT_DKD4(DRIFT_RI, KICK_RI, time_d, time_k, (DT));
	EVOLVE_SPLIT_CRI4(EVOLVE_RI, EVOLVE_CC, time_r, time_c, dt);
	//EVOLVE_SPLIT_CRI4(EVOLVE_CC, EVOLVE_RI, time_c, time_r, dt);

	clevel--;
	free_sys(c.next_cc);
}


#define CC2_4_SWITCH_THRESHOLD 10

void evolve_split_cc2_4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {
/*
 * Hybrid CC2_4 split.
 * Evolve system with 2nd order split.
 * If systems size is small enough in terms of particle numbers, evolve at 4th order
 */
  if (s.n > CC2_4_SWITCH_THRESHOLD) {
    struct sys c = zerosys, r = zerosys;
    clevel++;
    if (etime <= stime ||  dt==0 || clevel>=MAXLEVEL) endrun((char *)"timestep too small");
    split_cc(s, &c, &r, dt);
    //if (s.n != c.n) LOG_CC_SPLIT(&c, &r); // print out non-trivial splits
    if (IS_ZEROSYSs(c)) {
      deepsteps++;
      simtime+=dt;
    }
    // evolve all fast components, 1st time
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      evolve_split_cc2_4(*ci, stime, stime+dt/2, dt/2);
    }
    drift(r, stime+dt/2, dt/2); // drift r, 1st time
    // kick ci <-> cj
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      for (struct sys *cj = &c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
        if (ci != cj) {
          kick(*ci, *cj, dt);
          //kick(*cj, *ci, dt);
        }
      }
    }
    // kick c <-> rest
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      kick(r, *ci, dt);
      kick(*ci, r, dt);
    }
    kick(r, r, dt); // kick rest
    drift(r, etime, dt/2);  // drift r, 2nd time
    // evolve all fast components, 2nd time
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      evolve_split_cc2_4(*ci, stime+dt/2, etime, dt/2);
    }
    clevel--;
    free_sys(c.next_cc);
  } else {
    // try shared4 instead?
    // think of a heuristic for changing the time step parameter
    evolve_split_cc4(s, stime, etime, dt);
  }
}
