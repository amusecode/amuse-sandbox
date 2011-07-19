// Object containing (gravitational) force calculation

class FORCE
{
  mpreal softening;

  public:

  // constructors
  FORCE();
  FORCE(mpreal soft);

  // forces
  void gravity(STAR* s1, STAR* s2);
};

//////////////////////////////////////////
// Constructors
//////////////////////////////////////////
FORCE::FORCE()
{
  softening = "0";
}
FORCE::FORCE(mpreal soft)
{
  softening = soft;
}
//////////////////////////////////////////
// Forces
//////////////////////////////////////////
void FORCE::gravity(STAR* s1, STAR* s2)
{
  mpreal dx = s2->x - s1->x;
  mpreal dy = s2->y - s1->y;
  mpreal dz = s2->z - s1->z;

  mpreal dr2 = dx*dx + dy*dy + dz*dz + softening*softening;
  mpreal dr1 = sqrt(dr2);  
  mpreal dr3 = dr2*dr1;

  s1->add_ax( s2->m/dr3*dx );
  s1->add_ay( s2->m/dr3*dy );
  s1->add_az( s2->m/dr3*dz );
  s2->add_ax( s1->m/dr3*-dx );
  s2->add_ay( s1->m/dr3*-dy );
  s2->add_az( s1->m/dr3*-dz );
}

