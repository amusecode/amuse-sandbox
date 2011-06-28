#include "twobody/kepler.c"

void evolve_twobody(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {
  clevel++;
  if (etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small\n");
  if (s.n != 2) ENDRUN("two-body solver was called with sys.n=%u\n", s.n);
  // translate coordinates original frame to 2-body frame
  DOUBLE dpos[3];
  DOUBLE dvel[3];
  DOUBLE pos_cm[3];
  DOUBLE vel_cm[3];
  DOUBLE m1 = s.part->mass;
  DOUBLE m2 = s.last->mass;
  DOUBLE mtot = s.part->mass + s.last->mass;
  DOUBLE f1 = m1 / mtot;
  DOUBLE f2 = m2 / mtot;
  dpos[0] = s.part->pos[0] - s.last->pos[0];
  dpos[1] = s.part->pos[1] - s.last->pos[1];
  dpos[2] = s.part->pos[2] - s.last->pos[2];
  dvel[0] = s.part->vel[0] - s.last->vel[0];
  dvel[1] = s.part->vel[1] - s.last->vel[1];
  dvel[2] = s.part->vel[2] - s.last->vel[2];
  pos_cm[0] = (m1 * s.part->pos[0] + m2 * s.last->pos[0]) / mtot;
  pos_cm[1] = (m1 * s.part->pos[1] + m2 * s.last->pos[1]) / mtot;
  pos_cm[2] = (m1 * s.part->pos[2] + m2 * s.last->pos[2]) / mtot;
  vel_cm[0] = (m1 * s.part->vel[0] + m2 * s.last->vel[0]) / mtot;
  vel_cm[1] = (m1 * s.part->vel[1] + m2 * s.last->vel[1]) / mtot;
  vel_cm[2] = (m1 * s.part->vel[2] + m2 * s.last->vel[2]) / mtot;
  // call kepler solver
  kepler_set_position(dpos);
  kepler_set_velocity(dvel);
  set_mu(mtot);
  if (kepler_evolve(dt) != 0) ENDRUN("kepler_evolve() failed\n");
  kepler_get_position(dpos);
  kepler_get_velocity(dvel);
  // evolve center of mass for dt
  pos_cm[0] += vel_cm[0] * dt;
  pos_cm[1] += vel_cm[1] * dt;
  pos_cm[2] += vel_cm[2] * dt;
  // translate coordinates from 2-body frame to original frame
  s.part->pos[0] = pos_cm[0] + f1 * dpos[0];
  s.part->pos[1] = pos_cm[1] + f1 * dpos[1];
  s.part->pos[2] = pos_cm[2] + f1 * dpos[2];
  s.part->vel[0] = vel_cm[0] + f1 * dvel[0];
  s.part->vel[1] = vel_cm[1] + f1 * dvel[1];
  s.part->vel[2] = vel_cm[2] + f1 * dvel[2];
  s.last->pos[0] = pos_cm[0] - f2 * dpos[0];
  s.last->pos[1] = pos_cm[1] - f2 * dpos[1];
  s.last->pos[2] = pos_cm[2] - f2 * dpos[2];
  s.last->vel[0] = vel_cm[0] - f2 * dvel[0];
  s.last->vel[1] = vel_cm[1] - f2 * dvel[1];
  s.last->vel[2] = vel_cm[2] - f2 * dvel[2];
  // dummy test: evolve with shared global time step
  //evolve_unsplit4(s, stime, etime, dt, 1);
  // update statistics
  cecount[clevel]++;
  clevel--;
}
