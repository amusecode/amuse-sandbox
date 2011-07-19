// Object cluster containing stars

class CLUSTER
{
  mpreal t;
  int N;
  vector<STAR> star;

  DYNAMICS dynamics;
  FORCE force;
  mpreal a_max;

  public:

  // Constructors
  CLUSTER();
  CLUSTER(mpreal T, int n, vector<STAR> st);
  CLUSTER( string file );
  CLUSTER( string file, FORCE fo );

  // Set
  void set_t(mpreal T);

  // Get
  mpreal get_t();
  int get_N();
  mpreal get_a_max();
  STAR* get_pointer_to_star();

  // Calculate
  void calc_a();
  void leapfrog(mpreal dt);

  // Printers
  void print();
  void print( ofstream &data );
};

//////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////
CLUSTER::CLUSTER()
{
  t = "0";
  N = 0;
  star.clear();
  a_max = "0";
}
CLUSTER::CLUSTER(mpreal T, int n, vector<STAR> st)
{
  t = T;
  N = n;
  star = st;
  a_max = "0";
}
CLUSTER::CLUSTER( string file )
{
  star.clear();
  ifstream data;
  data.open( file.c_str() );
  if( !data )
  {
    cerr << "Could not open " << file << "!" << endl;
    exit(1);
  }
  else
  {
    data >> t >> N;
    mpreal m, x, y, z, vx, vy, vz;
    for(int i=0; i<N; i++)
    {
      data >> m >> x >> y >> z >> vx >> vy >> vz;
      STAR st(m, x, y, z, vx, vy, vz);
      star.push_back( st );
    }
  }
  data.close();
  a_max = "0";
}
CLUSTER::CLUSTER( string file, FORCE fo )
{
  star.clear();
  ifstream data;
  data.open( file.c_str() );
  if( !data )
  {
    cerr << "Could not open " << file << "!" << endl;
    exit(1);
  }
  else
  {
    data >> t >> N;
    mpreal m, x, y, z, vx, vy, vz;
    for(int i=0; i<N; i++)
    {
      data >> m >> x >> y >> z >> vx >> vy >> vz;
      STAR st(m, x, y, z, vx, vy, vz);
      star.push_back( st );
    }
  }
  data.close();
  a_max = "0";
  force = fo;
}
//////////////////////////////////////////////
// Set
//////////////////////////////////////////////
void CLUSTER::set_t(mpreal T)
{
  t = T;
}
//////////////////////////////////////////////
// Get
//////////////////////////////////////////////
mpreal CLUSTER::get_t()
{
  return t;
}
int CLUSTER::get_N()
{
  return N;
}
mpreal CLUSTER::get_a_max()
{
  return a_max;
}
STAR* CLUSTER::get_pointer_to_star()
{
  STAR *p = &star[0];
  return p;
}
//////////////////////////////////////////////
// Calculate
//////////////////////////////////////////////
void CLUSTER::calc_a()
{
  for(int i=0; i<N; i++)
  {
    star[i].reset_a();
  }
  for(int i=0; i<N-1; i++) 
  {
    for(int j=i+1; j<N; j++) 
    {
      force.gravity(&star[i], &star[j]);
    }
  }  
  a_max = "0";
  for(int i=0; i<N; i++) 
  {
    star[i].calc_a2mag();
    if( star[i].get_a2mag() > a_max )
    {
      a_max = star[i].get_a2mag();
    }
  }
  a_max = sqrt(a_max);
}
void CLUSTER::leapfrog(mpreal dt)
{
  for(int i=0; i<N; i++) dynamics.update_r(&star[i], dt);
  calc_a();
  for(int i=0; i<N; i++) dynamics.update_v(&star[i], dt);
}
//////////////////////////////////////////////
// Printers
//////////////////////////////////////////////
void CLUSTER::print()
{
  vector<STAR>::iterator st;
  cout << t << " " << N << endl;
  for(st=star.begin(); st!=star.end(); ++st)
  {
    cout << st->m << " " << st->x << " " << st->y << " " << st->z << " " << st->vx << " " << st->vy << " " << st->vz << endl;
  }
}
void CLUSTER::print( ofstream &data )
{
  vector<STAR>::iterator st;
  data << t << " " << N << endl;
  for(st = star.begin(); st != star.end(); ++st)
  {
    data << st->m << " " << st->x << " " << st->y << " " << st->z << " " << st->vx << " " << st->vy << " " << st->vz << endl;
  }
}

