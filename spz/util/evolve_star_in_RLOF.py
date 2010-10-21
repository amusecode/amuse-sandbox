# evolve a single star in isolation and with a limiting (Roche)
# radius.  using MESA. In the radius limiting case the evolution
# continues while keeping track of the mass loss, but the Roche radius
# is not self consistently adapted.
#    Leiden Observatory, 21 October 2010 
#    Nathan de Vries and Simon Portegies Zwart
#

import sys 
import numpy
from amuse.lab import *

#set up a single star of mass (mass) in MSun
#Returns a pointer to the star_particle in MESA
def setup_star_in_MESA(mass) :
    se = MESA()
    star_particle = Particle()
    star_particle.mass = mass
    se.initialize_module_with_default_parameters() 
    se.particles.add_particle(star_particle)
    se.initialize_stars()
    single_star = se.particles[0]
    return single_star 

# determine how much mass sticks out of radius rmax
# and returns dm
def determine_RLOF_mass_excess(star, rmax) :
    #obtain stellar structure (radius and mass)
    rprof = star.get_radius_profile()
    mfprof = star.get_mass_profile() # really the mass fraction
    #identify the shell where the radius exceeds rmax
    ri = numpy.searchsorted(rprof, rmax)
    # for the innermost mass-losing shell calculate the fraction of mass lost
    # determine the volume of the star that sticks out rmax
    dr_shell = rprof[ri]**2-rprof[ri-1]**2
    dr_RL = rprof[ri]**2-rmax**2
    drf = dr_RL/dr_shell
    dmf_shell = mfprof[ri] * drf 
    dm_shell = dmf_shell * star.mass
    # calculate the mass fraction of the star with rstar>rmax
    dm = dm_shell + (mfprof[ri:] * star.mass).sum() #summ the mass outside rmax
    # check the mass in the shells that are completely lost
    mtot = 0
    for i in numpy.arange(ri, star.get_number_of_zones().number) :
        mtot += mfprof[i]
    fm = mtot*star.mass/(dm-dm_shell)
    if abs(fm-1)>1.e-10 :
        print "mtot check=", fm
    return dm

# Apply Roche-lobe overflow
def RLOF(star, RRoche) :
    dm = determine_RLOF_mass_excess(star, RRoche)
    mass_scaling = (star.mass - dm) / star.mass
    star.mass -= dm

    # rescale stellar radius
    rprof = star.get_radius_profile()
    radial_scaling = RRoche/rprof[-1]
    rprof *= radial_scaling
    star.set_radius_profile(rprof)
    
    # rescale densities
    rho_prof = star.get_density_profile()
    rho_prof *= mass_scaling / radial_scaling**3
    star.set_density_profile(rho_prof)

    return dm
    
def evolve_star_in_RLOF(star, RRoche) :
    # evolve star until it overfills its Roche lobe
    while star.radius < RRoche :
        t_old = star.age
        star.evolve_one_step()
        print "Star: ", star

    # continue the evolution of the star through RLOF
    while star.radius > RRoche :
        dt = star.age - t_old
        dm = RLOF(star, RRoche)
        # dmdt = dm/dt
        # print "dmdt=", dm, dt, (dmdt).as_quantity_in(units.MSun/units.yr)
        t_old = star.age
        star.evolve_one_step()
        print star

# Evolve a single star precisely up to a certain age.
def evolve_star_to_time(star, t_end) :
    while star.age+star.time_step < t_end :
        star.evolve_one_step()
    star.time_step = t_end - star.age
    star.evolve_one_step()
    
if __name__ == '__main__' :

    t_end = -4412435911.67 | units.yr
    if len(sys.argv)>1 :
        m_init = float(sys.argv[1]) | units.MSun
        RRoche = float(sys.argv[2]) | units.RSun
        if len(sys.argv)>3 :
            t_end = float(sys.argv[3]) | units.Myr

    star = setup_star_in_MESA(m_init)
    if t_end.number>0  :
        evolve_star_to_time(star, t_end)
    else :
        evolve_star_in_RLOF(star, RRoche)
    print "star=", star

        


