// Object Star

class STAR
{
  public:

  // Variables
  mpreal m, x, y, z, vx, vy, vz;
  mpreal ax, ay, az, ax0, ay0, az0;
  mpreal r2_mag, v2_mag, a2_mag;

  // Constructors
  STAR();
  STAR(mpreal M, mpreal X, mpreal Y, mpreal Z, mpreal VX, mpreal VY, mpreal VZ);

  // Reset
  void reset_a();

  // Get
  mpreal get_a2mag();

  // Add
  void add_x(mpreal X);
  void add_y(mpreal Y);
  void add_z(mpreal Z);
  void add_vx(mpreal VX);
  void add_vy(mpreal VY);
  void add_vz(mpreal VZ);
  void add_ax(mpreal AX);
  void add_ay(mpreal AY);
  void add_az(mpreal AZ);

  // Calculate
  void calc_r2mag();
  void calc_v2mag();
  void calc_a2mag();
};

//////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////
STAR::STAR()
{
  m = "0";
  x = "0";
  y = "0";
  z = "0";
  vx = "0";
  vy = "0";
  vz = "0";
  ax = "0";
  ay = "0";
  az = "0";
  ax0 = "0";
  ay0 = "0";
  az0 = "0";
  r2_mag = "0";
  v2_mag= "0";
  a2_mag = "0";
} 
STAR::STAR(mpreal M, mpreal X, mpreal Y, mpreal Z, mpreal VX, mpreal VY, mpreal VZ)
{
  m = M;
  x = X;
  y = Y;
  z = Z;
  vx = VX;
  vy = VY;
  vz = VZ;
  ax = "0";
  ay = "0";
  az = "0";
  ax0 = "0";
  ay0 = "0";
  az0 = "0";
  r2_mag = "0";
  v2_mag = "0";
  a2_mag = "0";
}
//////////////////////////////////////////////
// Reset
//////////////////////////////////////////////
void STAR::reset_a()
{
  ax0 = ax;
  ay0 = ay;
  az0 = az;
  ax = "0";
  ay = "0";
  az = "0";
}
//////////////////////////////////////////////
// Get
//////////////////////////////////////////////
mpreal STAR::get_a2mag()
{
  return a2_mag;
}
//////////////////////////////////////////////
// Add
//////////////////////////////////////////////
void STAR::add_x(mpreal X)
{
  x += X;
}
void STAR::add_y(mpreal Y)
{
  y += Y;
}
void STAR::add_z(mpreal Z)
{
  z += Z;
}
void STAR::add_vx(mpreal VX)
{
  vx += VX;
}
void STAR::add_vy(mpreal VY)
{
  vy += VY;
}
void STAR::add_vz(mpreal VZ)
{
  vz += VZ;
}
void STAR::add_ax(mpreal AX)
{
  ax += AX;
}
void STAR::add_ay(mpreal AY)
{
  ay += AY;
}
void STAR::add_az(mpreal AZ)
{
  az += AZ;
}
//////////////////////////////////////////////
// Calc
//////////////////////////////////////////////
void STAR::calc_r2mag()
{
  r2_mag = x*x+y*y+z*z;
}
void STAR::calc_v2mag()
{
  v2_mag = vx*vx+vy*vy+vz*vz;
}
void STAR::calc_a2mag()
{
  a2_mag = ax*ax+ay*ay+az*az;
}

