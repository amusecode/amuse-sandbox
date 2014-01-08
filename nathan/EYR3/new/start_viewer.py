import sys
import os.path
import select
import numpy
import cPickle
from amuse.units import units, nbody_system
from amuse.io import read_set_from_file
from amuse.ext.particles_with_color import new_particles_with_blackbody_color
from amuse.community.asterisk.interface import Asterisk

def new_visualization(stars, gas):
    #creating visualization code
    converter = nbody_system.nbody_to_si(10.0 | units.parsec, stars.total_mass() + gas.total_mass())
    visualization = Asterisk(converter, redirection="none")

    #optional: change OpenGL perspective settings
    #visualization.set_field_of_view(45.0)
    #visualization.set_z_near(0.1 | units.parsec)
    #visualization.set_z_far(3000.0 | units.parsec)

    #initialize code (creates gui windows)
    visualization.initialize_code()
    visualization.parameters.use_octree_for_gas = True

    #optional: set the zoom and rotation of the visualization
    #visualization.parameters.rotation = (15, -15, 45)
    #visualization.parameters.camera_distance = 100 | units.parsec
    
    #add (now colored) particles to visualization
    visualization.star_particles.add_particles(stars)
    visualization.gas_particles.add_particles(gas)
    return visualization

def next_filename():
    next_filename.counter += 1
    return os.path.join("snapshots", "cluster_snapshot_{0:=06}_".format(next_filename.counter))
next_filename.counter = -1

def radius_from_luminosity(luminosity):
    return 0.02 * numpy.log(luminosity.value_in(units.LSun)) | units.parsec

if __name__ == "__main__":
    filename = next_filename()
    stars = read_set_from_file(filename+"stars.amuse", 'amuse')
    stars = new_particles_with_blackbody_color(stars)
    stars.radius = radius_from_luminosity(stars.luminosity)
    
    gas = read_set_from_file(filename+"gas.amuse", 'amuse')
    gas = new_particles_with_blackbody_color(gas)
    gas.alpha = 0.00001
    
    visualization = new_visualization(stars, gas)
    visualization.store_view(0|units.Myr)
    from_stars_to_viz = stars.new_channel_to(visualization.star_particles)
    from_gas_to_viz = gas.new_channel_to(visualization.gas_particles)
#~    print visualization.gas_particles.alpha
    
    filename = next_filename()
    while True:
        if not os.path.exists(filename+"gas.amuse"):
            print "Press 'Enter' to quit"
            input, o, e = select.select([sys.stdin], [], [], 3)
            if input:
                break
            continue
        
        next = read_set_from_file(filename+"stars.amuse", 'amuse', close_file=True)
        next.synchronize_to(stars)
        next.copy_values_of_attributes_to(["x", "y", "z", "temperature", "luminosity"], stars)
        stars.radius = radius_from_luminosity(stars.luminosity)
        stars.synchronize_to(visualization.star_particles)
        from_stars_to_viz.copy_attributes(["x", "y", "z", "radius", "red", "green", "blue"])
        
        next = read_set_from_file(filename+"gas.amuse", 'amuse', close_file=True)
        next.synchronize_to(gas)
        next.copy_values_of_attributes_to(["x", "y", "z", "u", "radius"], gas)
        gas.synchronize_to(visualization.gas_particles)
#~        from_gas_to_viz.copy_attributes(["x", "y", "z", "h_smooth", "red", "green", "blue"], ["x", "y", "z", "radius", "red", "green", "blue"])
        from_gas_to_viz.copy_attributes(["x", "y", "z", "radius", "red", "green", "blue"])
        
        with open(filename+"info.pkl", "rb") as infile:
            (tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, 
                total_feedback_energy, current_time
            ) = cPickle.load(infile)
        visualization.store_view("t={0}, feedback energy released: {1}".format(current_time, total_feedback_energy.as_quantity_in(units.erg)))
        filename = next_filename()
    visualization.stop()

