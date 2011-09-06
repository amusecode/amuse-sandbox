import sys

from amuse.community.simplex.interface import SimpleXInterface as mpi_interface

from amuse import datamodel
def rad_from_file(input_file):

  ids=[]
  rad = mpi_interface(number_of_workers = 1,redirection="none")
  rad.initialize_code()
  f = open(input_file,'r')
  lines = f.readlines()
  lines.pop(0)
  pid,x,y,z,nh,flux,xion=[],[],[],[],[],[],[]
  for line in lines:
    l = line.strip().split()
    if len(l) >= 7:
      pid.append(int(l[0]))
      x.append(float(l[1]))
      y.append(float(l[2]))
      z.append(float(l[3]))
      nh.append(float(l[4]))
      flux.append(float(l[5]))
      xion.append(float(l[6]))
      ids.append(int(l[0]))
  ids,n=rad.new_particle(x,y,z,nh,flux,xion)     
  
  rad.commit_particles()
  return rad,ids
  
# time unit=Myr
# intensity unit=1.e48 photon/sec
# mass unit= 6.757e64 amu 
# length unit= 13.20 kpc (SizeBox)
# density unit=0.001 amu/cm^3


if (__name__ == "__main__"):
  rad,ids=rad_from_file('vertices_64.txt')
#  rad.start_viewer()
  rad.evolve(0.1)
#  print rad.remove_particle(0)
#  n=rad.add_particle(1001,0.5,0.5,0.8,1.,50,0.99999)     
#  print n
#  rad.evolve(10.,1)
  results=rad.get_state(ids)
  del rad

  execfile('simplexplotje.py')
