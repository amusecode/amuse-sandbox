/*
 * Optimal kick Hamiltonian split (OK-split)
 */

struct force {
    struct particle *parti;
    struct particle *partj;
    FLOAT timestep;
};

struct forces {
    UINT n;
    struct force *forc;
    struct force *last;
};

struct forces zeroforces = {0, NULL, NULL};

DOUBLE ok_timestep_ij_fw(struct particle *i, struct particle *j) {

	FLOAT timestep;
	FLOAT dx[3],dr3,dr2,dr,dv[3],dv2,mu,vdotdr2,tau,dtau;
    timestep=HUGE_VAL;

    dx[0]=i->pos[0]-j->pos[0];
    dx[1]=i->pos[1]-j->pos[1];
    dx[2]=i->pos[2]-j->pos[2];
    dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;

    if(dr2>0) {
    	dr=sqrt(dr2);
    	dr3=dr*dr2;
        dv[0]=i->vel[0] - j->vel[0];
        dv[1]=i->vel[1] - j->vel[1];
        dv[2]=i->vel[2] - j->vel[2];
        vdotdr2=(dv[0]*dx[0]+dv[1]*dx[1]+dv[2]*dx[2])/dr2;
        dv2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
        mu=i->mass+j->mass;

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
		ENDRUN("negative timestep!\n");
	}

    return timestep;
}

DOUBLE ok_timestep_ij_bw(struct particle * i, struct particle * j) {

	FLOAT timestep;
	FLOAT dx[3],dr3,dr2,dr,dv[3],dv2,mu,vdotdr2,tau,dtau;
    timestep=HUGE_VAL;

    dx[0]=i->pos[0]-j->pos[0];
    dx[1]=i->pos[1]-j->pos[1];
    dx[2]=i->pos[2]-j->pos[2];
    dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;

    if(dr2>0) {
    	dr=sqrt(dr2);
    	dr3=dr*dr2;
        dv[0]=i->vel[0] - j->vel[0];
        dv[1]=i->vel[1] - j->vel[1];
        dv[2]=i->vel[2] - j->vel[2];
        vdotdr2=(dv[0]*dx[0]+dv[1]*dx[1]+dv[2]*dx[2])/dr2;
        dv2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
        mu=i->mass+j->mass;

#ifdef RATIMESTEP
        tau=RARVRATIO*dt_param/M_SQRT2*sqrt(dr3/mu);
        dtau=3/2.*tau*vdotdr2;
    	if(dtau<-1.) dtau=-1.;
        tau/=(1+dtau/2);
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

    timestep = -timestep;

    if (timestep > 0) {
		ENDRUN("positive timestep!\n");
	}

    return timestep;
}

static void ok_timestep_cpu(struct forces f, DOUBLE dt) {
	if (dt > 0) {
		for (UINT i = 0; i < f.n; i++) {
			//if (f.forc[i].timestep != HUGE_VAL) ENDRUN("timestep??");
			f.forc[i].timestep = ok_timestep_ij_fw(f.forc[i].parti, f.forc[i].partj);
		}
	} else {
		for (UINT i = 0; i < f.n; i++) {
			//if (f.forc[i].timestep != HUGE_VAL) ENDRUN("timestep??");
			f.forc[i].timestep = ok_timestep_ij_bw(f.forc[i].parti, f.forc[i].partj);
		}
	}
	tstep[clevel]++;
	tcount[clevel] += f.n;
}

#define LOG_FORCES(F) { \
	for (UINT i = 0; i < (F).n; i++) { \
		printf("%u\t%u\t%f\n", (F).forc[i].parti->id, (F).forc[i].partj->id, (F).forc[i].timestep); \
	} \
};

/*
 * split_ok_forces: split forces into smaller than dt, faster than dt
 */
static void ok_split(FLOAT dt, struct forces f, struct forces *slow, struct forces *fast) {
	//LOG("dt=%lf f.n=%u\n", dt, f.n);

	UINT i = 0;
	struct force *left, *right;
	left = f.forc;
	right = f.last;

	while (1) {
		if (i >= f.n) ENDRUN("forces split error 1\n");
		i++;
		while ((left->timestep < dt) && (left<right)) left++;
		while ((right->timestep >= dt) && (left<right)) right--;
		if (left < right) {
			SWAP(*left, *right, struct force);
		} else
			break;
	}

	if (left->timestep < dt) left++;
	slow->n = f.last - left + 1;
	fast->n = left - f.forc;

	if (fast->n == 1) {
		fast->n = 0;
		slow->n = f.n;
	}

	if (slow->n > 0) {
		slow->forc = f.forc + fast->n;
		slow->last = f.last;//slow->part+slow->n-1;
	}

	if (fast->n > 0) {
		fast->forc = f.forc;
		fast->last = f.forc + fast->n - 1;
	}
	if (fast->n + slow->n != f.n)
		ENDRUN("forces split error 2: fast->n=%u slow->n=%u f.n=%u\n", fast->n, slow->n, f.n);
	//for (i = 0; i < f.n; i++) f.forc[i].level = clevel;
}

struct forces ok_main_forces;

void ok_evolve_init(struct sys s) {
	UINT n_forces = s.n * s.n - s.n;
	ok_main_forces.forc = (struct force *) malloc(n_forces * sizeof(struct force));
	ok_main_forces.last = &(ok_main_forces.forc[n_forces - 1]);
	ok_main_forces.n = n_forces;

	// initialize pointers of the forces structure
	UINT k = 0;
	for (UINT i = 0; i < s.n; i++) {
		for (UINT j = 0; j < s.n; j++) {
			if (i != j) {
				ok_main_forces.forc[k].parti = &( s.part[i] );
				ok_main_forces.forc[k].partj = &( s.part[j] );
				k++;
			}
		}
	}
}

void ok_evolve_stop() {
	free(ok_main_forces.forc);
}

static void ok_kick(struct forces f, DOUBLE dt) {
	FLOAT dx[3],dr3,dr2,dr,acci;
	FLOAT acc[3];

	for (UINT i = 0; i < f.n; i++) {
		acc[0] = 0.;
		acc[1] = 0.;
		acc[2] = 0.;

		dx[0] = f.forc[i].parti->pos[0] - f.forc[i].partj->pos[0];
		dx[1] = f.forc[i].parti->pos[1] - f.forc[i].partj->pos[1];
		dx[2] = f.forc[i].parti->pos[2] - f.forc[i].partj->pos[2];
		dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + eps2;

		if (dr2 > 0) {
			dr = sqrt(dr2);
			dr3 = dr*dr2;
			acci = f.forc[i].partj->mass / dr3;

			f.forc[i].parti->vel[0] -= dt * dx[0] * acci;
			f.forc[i].parti->vel[1] -= dt * dx[1] * acci;
			f.forc[i].parti->vel[2] -= dt * dx[2] * acci;
		}
	}

	kstep[clevel]++;
	kcount[clevel] += f.n;
}

#define IS_ZEROFORCES(F) (((F).n == 0) && ((F).forc == NULL) && ((F).last == NULL))

void ok_evolve2(struct sys s, struct forces f, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {

	if (IS_ZEROFORCES(f) && clevel == -1) { f = ok_main_forces; }
	clevel++;
    if ((etime <= stime) || (dt == 0) || (clevel >= MAXLEVEL))
      ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u", etime, stime, dt, clevel);

    // all particles are drifted together
    if (f.n == 0) {
    	deepsteps++;
    	simtime += dt;
    	drift(s, etime, dt);
    	clevel--;
    	return;
    }

    //printf("forces before timestep_cpu_ok():\n");
	//LOG_FORCES(f);

	if (calc_timestep) ok_timestep_cpu(f, dt);
	//LOG("forces f before split_ok_forces():\n");
	//LOG_FORCES(f);

	struct forces slowf = zeroforces, fastf = zeroforces;
	ok_split((FLOAT) dt, f, &slowf, &fastf);
	//LOG("forces f after split_ok_forces():\n");
	//LOG_FORCES(f);
	//LOG("fast after split_ok_forces():\n");
	//LOG_FORCES(fastf);
	//LOG("slow after split_ok_forces():\n");
	//LOG_FORCES(slowf);

	ok_evolve2(s, fastf, stime, stime+dt/2, dt/2, 0);
	ok_kick(slowf, dt);
	ok_evolve2(s, fastf, stime+dt/2, etime, dt/2, 1);

    clevel--;
}

void ok_evolve4(struct sys s, struct forces f, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {

	if (IS_ZEROFORCES(f) && clevel == -1) { f = ok_main_forces; }
	clevel++;
    if (((dt > 0) && (etime <= stime)) || ((dt < 0) && (etime >= stime)) || // generalization of: etime <= stime
    	(dt == 0) || (clevel >= MAXLEVEL))
    		ENDRUN("timestep too small: stime=%Lf etime=%Lf dt=%Lf clevel=%u\n", stime, etime, dt, clevel);

    // all particles are drifted together
    if (f.n == 0) {
    	deepsteps++;
    	simtime += dt;
    	drift(s, etime, dt);
    	clevel--;
    	return;
    }

	if (calc_timestep) ok_timestep_cpu(f, dt);
	//LOG("ok_evolve4() dt=%Lf\n", dt);
	//LOG_FORCES(f);

	struct forces slowf = zeroforces, fastf = zeroforces;
	if (dt > 0) {
		ok_split((FLOAT) dt, f, &slowf, &fastf);
	} else {
		ok_split((FLOAT) dt, f, &fastf, &slowf);
	}

	DOUBLE atime = stime, btime = stime;

	#define EVOLVE_A(TIME, DT) \
		ok_evolve4(s, fastf, (TIME), (TIME) + (DT), (DT), 1); \
		(TIME) += (DT);

	#define EVOLVE_B(TIME, DT) \
		ok_kick(slowf, (DT)); \
		(TIME) += (DT);

	// 2nd order split, forward time, tests operator split macros
	//EVOLVE_A(atime, 0.5*dt);
	//EVOLVE_B(btime, 1.0*dt);
	//EVOLVE_A(atime, 0.5*dt);

	// 2nd order split, backward time, tests operator split macros
	//EVOLVE_A(atime, -0.5*ABS(dt));
	//EVOLVE_B(btime, -1.0*ABS(dt));
	//EVOLVE_A(atime, -0.5*ABS(dt));

	// 4-th order symmetric operator split
	#define A1 0.0792036964311957
	#define A2 0.353172906049774
	#define A3 -.0420650803577195
	#define A4 (1 - 2*((A1) + (A2) + (A3))) // A4=0.21937695575349958
	#define B1 0.209515106613362
	#define B2 -.143851773179818
	#define B3 (0.5 - (B1) - (B2)) // B3=0.43433666656645598

	EVOLVE_A(atime, A1*dt);
	EVOLVE_B(btime, B1*dt);
	EVOLVE_A(atime, A2*dt);
	EVOLVE_B(btime, B2*dt);
	EVOLVE_A(atime, A3*dt);
	EVOLVE_B(btime, B3*dt);
	EVOLVE_A(atime, A4*dt);
	EVOLVE_B(btime, B3*dt);
	EVOLVE_A(atime, A3*dt);
	EVOLVE_B(btime, B2*dt);
	EVOLVE_A(atime, A2*dt);
	EVOLVE_B(btime, B1*dt);
	EVOLVE_A(atime, A1*dt);

    clevel--;
}

