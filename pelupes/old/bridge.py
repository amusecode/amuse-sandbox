# kick system sys1 according to gravity func get_gravity
def kick_system(sys1, get_gravity, dt,ids=None):
  if(ids is None):
    all_ids=[]
    error, current = sys1.get_index_of_first_particle()
    while error == 0:
      all_ids.append(current)
      error, current_index = instance.get_index_of_next_particle(current_index)
    ids=all_ids  
  
  for pid in ids:
    (mass,radius,x,y,z,vx,vy,vz),error=sys1.get_state(pid)
    ax,ay,az,error=get_gravity(radius,x,y,z)
    vx=vx+dt*ax
    vy=vy+dt*ay
    vz=vz+dt*az
    sys1.set_state(pid,mass,radius,x,y,z,vx,vy,vz)

def evolve_w_kick(tend,sys,perturber,ids=None):

  get_gravity=perturber.get_gravity_at_point
  tnow,err=sys.get_time()
  dt=tend-tnow
  perturber.evolve(tnow)
  kick_system(sys,get_gravity,dt/2,ids) 

  sys.evolve(tend)

  perturber.evolve(tend)
  kick_system(sys,get_gravity,dt/2,ids)  

