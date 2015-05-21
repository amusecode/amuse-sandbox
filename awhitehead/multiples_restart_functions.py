import os
import random
import pickle
import numpy

from amuse.couple import multiples
from amuse import datamodel
from amuse import io

#We save both the python set of stars and the set of stars as seen by gravity.
#Extra_stuff is whatever additional information the user wants,
#assuming it can be pickled without breaking...
#If backup == 1, the code will write a duplicate set of save files with .backup before other extensions.  The idea here is to
#remove the possibility of having only one set of write files get destroyed by a system failure during the writing steps.
def write_state_to_file(time, stars_python,gravity_code, multiples_code, write_file, extra_stuff=None,cp_hist=False, backup = 0 ):
    print "Writing state to write file: ", write_file,"\n"
    if write_file is not None:
        particles = gravity_code.particles.copy()
        write_channel = gravity_code.particles.new_channel_to(particles)
        write_channel.copy_attribute("index_in_code", "id")
        bookkeeping = {'neighbor_veto': multiples_code.neighbor_veto,
            'neighbor_distance_factor': multiples_code.neighbor_distance_factor,
            'multiples_external_tidal_correction': multiples_code.multiples_external_tidal_correction,
            'multiples_integration_energy_error': multiples_code.multiples_integration_energy_error,
            'multiples_internal_tidal_correction': multiples_code.multiples_internal_tidal_correction,
            'model_time': multiples_code.model_time,
            'root_index': multiples.root_index
        }

        for root, tree in multiples_code.root_to_tree.iteritems():
            try:
                #multiples.print_multiple_simple(tree,kep)
                root_in_particles = root.as_particle_in_set(particles)
                subset = tree.get_tree_subset().copy()
                root_in_particles.components = subset
            except AttributeError as e:
                try:
                    print "Attribute error on writing: %s" % e
                    multiples.print_multiple_simple(tree,kep)
                except:
                    pass

        io.write_set_to_file(particles,write_file+".stars.hdf5",'hdf5',version='2.0', append_to_file=False, copy_history=cp_hist)
        io.write_set_to_file(stars_python,write_file+".stars_python.hdf5",'hdf5',version='2.0', append_to_file=False, copy_history=cp_hist)
        config = {'time' : time,
            'py_seed': pickle.dumps(random.getstate()),
            'numpy_seed': pickle.dumps(numpy.random.get_state()),
        }
        with open(write_file + ".conf", "wb") as f:
            pickle.dump(config, f)

        with open(write_file + ".bookkeeping", "wb") as f:
            pickle.dump(bookkeeping, f)

        if extra_stuff is not None:
            with open(write_file + ".extra_stuff", "wb") as f:
                pickle.dump(extra_stuff, f)

        print "\nState successfully written to:  ", write_file
        print time
        if backup == 1:
            io.write_set_to_file(particles,write_file+".backup.stars.hdf5",'hdf5',version='2.0', append_to_file=False, copy_history=cp_hist)
            io.write_set_to_file(stars_python,write_file+".backup.stars_python.hdf5",'hdf5',version='2.0', append_to_file=False, copy_history=cp_hist)
            with open(write_file + ".backup.conf", "wb") as f:
                pickle.dump(config, f)
            with open(write_file + ".backup.bookkeeping", "wb") as f:
                pickle.dump(bookkeeping, f)
            if extra_stuff is not None:
                with open(write_file + ".backup.extra_stuff", "wb") as f:
                    pickle.dump(extra_stuff, f)
            print "\nBackup write completed.\n"
#restores the state of the system from the data files created by write_state_from_file.  Options loaded separately.
#initializes and returns a multiples code
def read_state_from_file(restart_file, gravity_code,resolve_collision_code_creation_function, kep):

    stars = io.read_set_from_file(restart_file+".stars.hdf5",'hdf5',version='2.0').copy()
    stars_python = io.read_set_from_file(restart_file+".stars_python.hdf5",'hdf5',version='2.0').copy()
    with open(restart_file + ".bookkeeping", "rb") as f:
        bookkeeping = pickle.load(f)
    if os.path.isfile(restart_file + ".extra_stuff"):
        with open(restart_file + ".extra_stuff", "rb") as f:
            extra_stuff = pickle.load(f)
    else:
        extra_stuff = None
    print bookkeeping
    root_to_tree = {}
    for root in stars:
        if hasattr(root, 'components') and not root.components is None:
            root_to_tree[root] = datamodel.trees.BinaryTreeOnParticle(root.components[0])
    gravity_code.set_begin_time(bookkeeping['model_time'])
    gravity_code.particles.add_particles(stars)
    gravity_code.commit_particles()


    multiples_code = multiples.Multiples(gravity_code, resolve_collision_code_creation_function, kep)
    #multiples_code.neighbor_distance_factor = 1.0
    #multiples_code.neighbor_veto = False
    #multiples_code.neighbor_distance_factor = 2.0
    #multiples_code.neighbor_veto = True
    multiples_code.neighbor_distance_factor = bookkeeping['neighbor_distance_factor']
    multiples_code.neighbor_veto = bookkeeping['neighbor_veto']
    multiples_code.multiples_external_tidal_correction = bookkeeping['multiples_external_tidal_correction']
    multiples_code.multiples_integration_energy_error = bookkeeping['multiples_integration_energy_error']
    multiples_code.multiples_internal_tidal_correction = bookkeeping['multiples_internal_tidal_correction']
    multiples.root_index = bookkeeping['root_index']
    multiples_code.root_to_tree = root_to_tree
    #multiples_code.set_model_time = bookkeeping['model_time']
    
    with open(restart_file + ".conf", "rb") as f:
        config = pickle.load(f)
    random.setstate(pickle.loads(config["py_seed"]))
    numpy.random.set_state(pickle.loads(config["numpy_seed"]))

    return stars_python, bookkeeping['model_time'], multiples_code, extra_stuff
