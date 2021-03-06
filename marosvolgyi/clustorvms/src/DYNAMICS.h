// Object containing leapfrog integrator

class DYNAMICS
{
  public:

  void update_r(STAR* s, mpreal dt);
  void update_v(STAR* s, mpreal dt);
};

void DYNAMICS::update_r(STAR* s, mpreal dt)
{
  s->add_x( s->vx*dt + "0.5"*s->ax*dt*dt );
  s->add_y( s->vy*dt + "0.5"*s->ay*dt*dt );
  s->add_z( s->vz*dt + "0.5"*s->az*dt*dt );
}
void DYNAMICS::update_v(STAR* s, mpreal dt)
{
  s->add_vx( "0.5"*(s->ax+s->ax0)*dt );
  s->add_vy( "0.5"*(s->ay+s->ay0)*dt );
  s->add_vz( "0.5"*(s->az+s->az0)*dt );
}

