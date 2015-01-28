######################################################

## This is a test to see whether Multiples can detect
## a new dynamically formed binary. 

## We start with the 3-body equilateral triangle
## configuration, which are 3 single stars. This 
## configuration falls apart and after a phase of 
## resonance, a binary and escaper are formed. 

## We output as a function of time, whether Multiples
## thinks there is a binary, and we output the result
## of a manual check for the presence of a binary.  

######################################################

import math

from amuse.datamodel import Particles
from amuse.units import nbody_system

from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler
from amuse.community.ph4.interface import ph4

from amuse.couple import encounters

######################################################

def new_smalln():
    code = SmallN()
    code.initialize_code()
    code.parameters.set_defaults()
    return code

def new_kepler():
    code = Kepler()
    code.initialize_code()
    code.parameters.set_defaults()
    return code

def new_ph4():
    code = ph4()
    code.initialize_code()
    code.parameters.set_defaults()
    return code

######################################################

def get_triangle():
    star = Particles(3)

    star[0].mass = 1. | nbody_system.mass
    star[0].x = 0. | nbody_system.length 
    star[0].y = 0. | nbody_system.length
    star[0].z = 0. | nbody_system.length
    star[0].vx = -0.5 | nbody_system.speed 
    star[0].vy = math.sqrt(1. - 0.5*0.5) | nbody_system.speed
    star[0].vz = 0. | nbody_system.speed

    star[1].mass = 1. | nbody_system.mass
    star[1].x = 1.0 | nbody_system.length  
    star[1].y = 0. | nbody_system.length
    star[1].z = 0. | nbody_system.length
    star[1].vx = -0.5 | nbody_system.speed  
    star[1].vy = -math.sqrt(1. - 0.5*0.5) | nbody_system.speed
    star[1].vz = 0. | nbody_system.speed 

    star[2].mass = 1. | nbody_system.mass
    star[2].x = 0.5 | nbody_system.length 
    star[2].y = math.sqrt(1. - 0.5*0.5) | nbody_system.length
    star[2].z = 0. | nbody_system.length
    star[2].vx = 1.0 | nbody_system.speed  
    star[2].vy = 0. | nbody_system.speed
    star[2].vz = 0. | nbody_system.speed   
 
    star[0].radius = 0. | nbody_system.length
    star[1].radius = 0. | nbody_system.length
    star[2].radius = 0. | nbody_system.length

    return star

######################################################

def check_number_of_binaries(N, m, pos, vel):
    Nbin = 0

    i=0
    while i<N-1:
        j=i+1
        while j<N:

            mu  = (m[i]+m[j]).number
            dx  = (pos[i][0]-pos[j][0]).number
            dy  = (pos[i][1]-pos[j][1]).number
            dz  = (pos[i][2]-pos[j][2]).number
            dvx = (vel[i][0]-vel[j][0]).number
            dvy = (vel[i][1]-vel[j][1]).number
            dvz = (vel[i][2]-vel[j][2]).number

            dr2 = dx**2+dy**2+dz**2
            dr = math.sqrt(dr2)

            dv2 = dvx**2+dvy**2+dvz**2

            # Negative energy
            e = 0.5*dv2 - mu/dr
            if e < 0.:

                lx = dy*dvz - dz*dvy
                ly = dz*dvx - dx*dvz
                lz = dx*dvy - dy*dvx
                l2 = lx**2+ly**2+lz**2

                D = mu**2+2.*e*l2
                if D >= 0.:
                    rp = (-mu + math.sqrt(D))/(2.*e)
                    ra = (-mu - math.sqrt(D))/(2.*e) 
                    a = 0.5*(ra+rp)

                    # If isolated and nearest neighbours
                    if i == 0 and j == 1:
                        index_nn = 2
                    elif i == 0 and j == 2:
                        index_nn = 1
                    else:
                        index_nn = 0 

                    dxi  = (pos[i][0]-pos[index_nn][0]).number
                    dyi  = (pos[i][1]-pos[index_nn][1]).number
                    dzi  = (pos[i][2]-pos[index_nn][2]).number                    
                    dr2i = dxi**2+dyi**2+dzi**2
                    dri = math.sqrt(dr2i)

                    dxj  = (pos[j][0]-pos[index_nn][0]).number
                    dyj  = (pos[j][1]-pos[index_nn][1]).number
                    dzj  = (pos[j][2]-pos[index_nn][2]).number                    
                    dr2j = dxj**2+dyj**2+dzj**2
                    drj = math.sqrt(dr2j)

                    if dri/a > 3. and drj/a > 3.:
                        Nbin += 1  

            j += 1
        i += 1

    return Nbin

######################################################

def run():
    # Generate initial conditions
    stars = get_triangle()

    # Setup dynamics 
    grav = new_ph4()

    encounter_code = encounters.HandleEncounter(
        kepler_code =  new_kepler(),
        resolve_collision_code = new_smalln(),
        interaction_over_code = None,
    )
    multiples_code = encounters.Multiples(
        gravity_code = grav,
        handle_encounter_code = encounter_code,
    )

    multiples_code.particles.add_particles(stars)
    multiples_code.commit_particles()

    # Setup simulation parameters and diagnostics
    t = 0. | nbody_system.time
    t_end = 100. | nbody_system.time
    dt = 1.0 | nbody_system.time

    stars = multiples_code.all_singles
    N = len(stars)
    m = stars.mass
    pos = stars.position
    vel = stars.velocity

    print 't=', t.number, 'N=', len(multiples_code.particles), 'Ns=', len(multiples_code.singles), 'Nb=', len(multiples_code.binaries), 'Nsib=', len(multiples_code.singles_in_binaries), 'Nm=', len(multiples_code.multiples), 'Ncom=', len(multiples_code.components_of_multiples), 'Nbin_manual=', check_number_of_binaries(N, m, pos, vel)

    # Start the integration
    while t<t_end:
        t += dt

        multiples_code.evolve_model(t)

        stars = multiples_code.all_singles
        N = len(stars)
        m = stars.mass
        pos = stars.position
        vel = stars.velocity

        print 't=', t.number, 'N=', len(multiples_code.particles), 'Ns=', len(multiples_code.singles), 'Nb=', len(multiples_code.binaries), 'Nsib=', len(multiples_code.singles_in_binaries), 'Nm=', len(multiples_code.multiples), 'Ncom=', len(multiples_code.components_of_multiples), 'Nbin_manual=', check_number_of_binaries(N, m, pos, vel)

    # Cleanup
    grav.cleanup_code()
    grav.stop()  

######################################################

if __name__ == "__main__":

    run()


