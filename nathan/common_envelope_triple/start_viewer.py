import sys
import os.path
import select
import numpy
from amuse.units import units, nbody_system
from amuse.io import read_set_from_file
from amuse.ext.particles_with_color import new_particles_with_blackbody_color
from amuse.community.asterisk.interface import Asterisk

def new_visualization(core, gas):
    converter = nbody_system.nbody_to_si(10.0 | units.AU, core.total_mass() + gas.total_mass())
    visualization = Asterisk(converter)
    visualization.initialize_code()
    #optional: set the zoom and rotation of the visualization
    #visualization.parameters.rotation = (15, -15, 45)
    #visualization.parameters.camera_distance = 100 | units.parsec
    visualization.marker_particles.add_particles(core)
    visualization.gas_particles.add_particles(gas)
    return visualization

def next_filename():
    next_filename.counter += 1
    return os.path.join("snapshots", "hydro_giant_{0:=04}_".format(next_filename.counter))
next_filename.counter = -1

if __name__ == "__main__":
    filename = next_filename()
    core = read_set_from_file(filename+"dm.amuse", 'amuse')
    core.red, core.green, core.blue = [1.0, 0.0, 0.0]
    
    gas = read_set_from_file(filename+"gas.amuse", 'amuse', close_file=True)[-9500:].copy()
    gas = new_particles_with_blackbody_color(gas)
    
    visualization = new_visualization(core, gas)
    visualization.store_view("Main sequence stars only")
    from_gas_to_viz = gas.new_channel_to(visualization.gas_particles)
    
    filename = next_filename()
    while True:
        if not os.path.exists(filename+"gas.amuse"):
            print "Press 'Enter' to quit"
            input, o, e = select.select([sys.stdin], [], [], 3)
            if input:
                break
            continue
        
        next = read_set_from_file(filename+"dm.amuse", 'amuse', close_file=True)
        next.copy_values_of_attributes_to(["x", "y", "z"], visualization.marker_particles)
        
        next = read_set_from_file(filename+"gas.amuse", 'amuse', close_file=True)#[-9500:]
        next.synchronize_to(gas)
        next.copy_values_of_attributes_to(["x", "y", "z", "u", "radius"], gas)
        gas.synchronize_to(visualization.gas_particles)
        from_gas_to_viz.copy_attributes(["x", "y", "z", "radius", "red", "green", "blue"])
        
        visualization.store_view()
        filename = next_filename()
    visualization.stop()

