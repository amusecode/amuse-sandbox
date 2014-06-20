from amuse.lab import *

def make_triple(inner_primary, inner_secondary, outer_star, inner_semimajor_axis, outer_semimajor_axis, inner_eccentricity=0.0, outer_eccentricity=0.0):

    stars = Particles(3)
    stars.is_star = True
    stars.add_particle(inner_primary)
    stars.add_particle(inner_secondary)
    stars.add_particle(outer_star)
    # this goes wrong, there are now 6 particles
    
    bins = Particles(2)
    bins.is_star = False
    bins[0].child1 = stars[0]
    bins[0].child2 = stars[1]
    bins[0].semimajor_axis = inner_semimajor_axis
    bins[1].semimajor_axis = outer_semimajor_axis
    bins[0].eccentricity = inner_eccentricity
    bins[1].eccentricity = outer_eccentricity
    bins[1].child1 = stars[2]
    bins[1].child2 = bins[0]
    

    triples = Particles(1)
    triples[0].inner_binary = bins[0]
    triples[0].outer_binary = bins[1]
    for x in bins:
        x.mass = x.child1.mass + x.child2.mass

    return triples

def evolve_binary(binary):
    print "evolve binary"
    return
def evolve_center_of_mass(binary):
    print "evolve center of mass"
    return

def evolve_triple(triple, time):
    
    print 'inner binary'
    if is_binary(triple.inner_binary):
        evolve_binary(triple.inner_binary)
    else:
        evolve_center_of_mass(triple.inner_binary)

    print 'outer binary'
    if is_binary(triple.outer_binary):
        evolve_binary(triple.inner_binary)
    else:
        evolve_center_of_mass(triple.inner_binary)



if __name__ == '__main__':
    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.yr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    p = Particles(mass=[1,2,3]|units.MSun)
    time = 100|units.yr
    
    
    triples = make_triple(p[0], p[1], p[2], 1|units.AU, 10|units.AU)
    print "Triples:", triples
    triples[0].inner_binary.semimajor_axis += 10|units.RSun
    triples[0].inner_binary.eccentricity = 0.6
    #write_set_to_file(coms, "file.data", "hdf5", version=2)
    print "Triples:", triples.inner_binary
    print "Triples:", triples.outer_binary
    
    print triples[0].outer_binary.child1.is_star
    print triples[0].outer_binary.child2.is_star
    
    