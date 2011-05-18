DOUBLE system_kinetic_energy(struct sys s) {
  DOUBLE e = 0.;
  for (UINT i = 0; i < s.n; i++)
    e += 0.5 * s.part[i].mass * (
        s.part[i].vel[0] * s.part[i].vel[0] +
        s.part[i].vel[1] * s.part[i].vel[1] +
        s.part[i].vel[2] * s.part[i].vel[2]);
  return e;
}

DOUBLE system_potential_energy(struct sys s) {
  DOUBLE e = 0.;
  for (UINT i = 0; i < s.n; i++)
    e += s.part[i].mass * s.part[i].pot;
  return e / 2;
}

DOUBLE system_linear_momentum_x(struct sys s) {
  DOUBLE p = 0.;
  for (UINT i = 0; i < s.n; i++) {
    p += s.part[i].mass * s.part[i].vel[0];
    LOG("px=%lE\n", p);
  }
  return p;
}

DOUBLE system_angular_momentum_x(struct sys s) {
  UINT i;
  DOUBLE l = 0.;
  for (i = 0; i < s.n; i++) {
    l += s.part[i].mass * (s.part[i].pos[1]*s.part[i].vel[2] + s.part[i].pos[2]*s.part[i].vel[1]);
    LOG("Lx=%lE\n", l);
  }
  return l;
}

DOUBLE system_center_of_mass_x(struct sys s) {
  DOUBLE x = 0.;
  for (UINT i = 0; i < s.n; i++)
    x += s.part[i].mass * s.part[i].pos[0];
  return x;
}

