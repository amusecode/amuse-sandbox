"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
import sys
import numpy
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from time import sleep

def make_map(disk, N, L):

    disk.h_smooth = 0 | units.RSun
    disk.u = 0 | (nbody_system.length/nbody_system.time)**2
#    Rmax = disk.position.lengths().amax()
    Rmax = L
    converter=nbody_system.nbody_to_si(disk.mass.sum(), Rmax)
    hydro=Fi(converter)

    hydro.parameters.use_hydro_flag=True
    hydro.parameters.radiation_flag=False
    hydro.parameters.self_gravity_flag=True
    hydro.parameters.gamma=1.
    hydro.parameters.isothermal_flag=True
    hydro.parameters.integrate_entropy_flag=False
    hydro.parameters.timestep=0.0001| units.day
    hydro.parameters.periodic_box_size=2 | units.AU

    hydro.gas_particles.add_particles(disk)
    grid = Grid.create((N,N,1), [2*Rmax, 2*Rmax, 2*Rmax])

    grid.position -= Rmax*[1,1,1]
    grid.vx = 0 | units.kms
    grid.vy = 0 | units.kms
    grid.vz = 0 | units.kms
    
    rho,rhovx,rhovy,rhovz,rhoe=hydro.get_hydro_state_at_point(grid.x.flatten(),grid.y.flatten(),grid.z.flatten(),grid.vx.flatten(),grid.vy.flatten(),grid.vz.flatten())
    rho=rho.reshape((N,N))

    hydro.stop()
    return numpy.transpose(rho)


def main(filename = "nbody.hdf5", lim=-1):
    L = lim 
#    pyplot.ion()
    storage = store.StoreHDF(filename,"r")
    stars = storage.load()
#    lim = max(stars.x).value_in(stars.x.unit)
    m =  0.05+10.0*stars.mass/max(stars.mass)
    i = 0
    time = 0 | units.day
    dt = 1./240. | units.day
    for si in stars.history:
        i+=1
        time += dt
        si = si.copy()
        si.move_to_center()
        pyplot.figure(figsize=(8,8))
        rho = make_map(si[2:], 100, L| units.RSun)

        pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
        extent=[-L,L,-L,L]) #,vmin=10,vmax=15)    

        pyplot.title("Accretor at t="+str(time.value_in(units.day))+"days")
#    pyplot.savefig('test.png')

#        pyplot.title("Cluster at t="+str(time))
        print "time = ", time+dt
#        scatter(si[0].vx, si[0].vy, s=100, c='r')
#        scatter(si[1].vx, si[1].vy, s=100, c='r')
#        scatter(si.vx, si.vy, s=m)
        pyplot.xlabel("X [$R_\odot$]")
        pyplot.ylabel("Y [$R_\odot$]")
#        if lim>0:
#            pyplot.xlim(-lim, lim)
#            pyplot.ylim(-lim, lim)
#        lim = 3.e+20
#        pyplot.xlim(-lim, lim)
#        pyplot.ylim(-lim, lim)
#        pyplot.draw()
        pyplot.savefig('map/map_'+'%6.6i.png'%i)
#        pyplot.savefig(filename)
#        sleep(10)
#        pyplot.cla()
#    pyplot.show()
#    sleep(100)

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "nbody.hdf5",
                      help="output filename [nbody.hdf5]")
    result.add_option("-l", dest="lim", type="float", default = -1,
                      help="boxsize")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


