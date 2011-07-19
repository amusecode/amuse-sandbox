// Clock object which handles timers and timestep size

class CLOCK
{
  mpreal t, t_end;
  mpreal dt, dt_min, dt_max, dv_max;
  mpreal t_print, dt_print, t_print0;

  mpreal t_cpu, t_progress;
  struct timeval Tvalue;
  struct timezone dummy;
  bool timerStarted;

  public:

  // constructors
  CLOCK();
  CLOCK(mpreal T_BEGIN, mpreal T_SIM, mpreal DT_MIN, mpreal DT_MAX, mpreal DV_MAX, mpreal DT_PRINT);

  // get
  mpreal get_t();
  mpreal get_dt();
  mpreal get_t_print();

  // runtime clock
  void start_timer();
  void stop_timer();
  mpreal get_timer();
  mpreal get_progress();

  // nbody clock
  int alarm();
  int to_print();
  void abort();
  void tick();
  void calc_dt( mpreal a_max );
};

///////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////
CLOCK::CLOCK()
{
  t = "0";
  t_end = "1";
  dt_min = "1e-4";
  dt_max = "1e-1";
  dv_max = "1e-2";
  dt_print = "1e-2";
  t_print = "1e-2";
  t_cpu = "0";
  dt = dt_min;
  t_progress = "0";
}
CLOCK::CLOCK(mpreal T_BEGIN, mpreal T_SIM, mpreal DT_MIN, mpreal DT_MAX, mpreal DV_MAX, mpreal DT_PRINT)
{
  t = T_BEGIN;
  t_end = T_BEGIN+T_SIM;
  dt_min = DT_MIN;
  dt_max = DT_MAX;
  dv_max = DV_MAX;
  dt_print = DT_PRINT;
  t_print = T_BEGIN + DT_PRINT;
  t_cpu = "0";
  dt = dt_min;
  t_progress = "0";
}
///////////////////////////////////////////////////////////////
// Get
///////////////////////////////////////////////////////////////
mpreal CLOCK::get_t()
{
  return t;
}
mpreal CLOCK::get_dt()
{
  return dt;
}
mpreal CLOCK::get_t_print()
{
  return t_print0;
}
///////////////////////////////////////////////////////////////
// runtime clock
///////////////////////////////////////////////////////////////
void CLOCK::start_timer()
{
  gettimeofday(&Tvalue,&dummy);
  timerStarted = true;
}
void CLOCK::stop_timer()
{
  if(timerStarted == false) t_cpu = "0";
  struct timeval Tvalue2;
  struct timezone dummy2;
  gettimeofday(&Tvalue2,&dummy2);
  mpreal startTime =  ((mpreal) Tvalue.tv_sec +"1.e-6"*((mpreal) Tvalue.tv_usec));
  mpreal endTime =  ((mpreal) Tvalue2.tv_sec +"1.e-6"*((mpreal) Tvalue2.tv_usec));
  timerStarted = false;
  t_cpu = endTime-startTime;      
} 
mpreal CLOCK::get_timer()
{
  return t_cpu;
}
mpreal CLOCK::get_progress()
{
  return t/t_end*"100.0";
}
///////////////////////////////////////////////////////////////
// nbody clock
///////////////////////////////////////////////////////////////
int CLOCK::alarm()
{
  if(t_print <= t_end+"0.5"*dt_print) return 0;
  else return 1;
}
int CLOCK::to_print()
{
  if(t > t_print)
  {
    t_print0 = t_print;
    t_print += dt_print;
    return 1;
  }
  else return 0;
}
void CLOCK::abort()
{
  t_print = t_end + dt_print;
}
void CLOCK::tick()
{
  t += dt;
}
void CLOCK::calc_dt( mpreal a_max )
{
  dt = dv_max/a_max;
  if(dt < dt_min) dt = dt_min;
  else if(dt > dt_max) dt = dt_max;
}

