import h5py
import os,sys
from amuse.lab import *
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4

snipfilename = sys.argv[1]

def readgadgethdfparticles(ggparticles, ptype="dm"):
    particles                 = Particles(len(ggparticles['ParticleIDs']))
    particles.position        = ggparticles['Coordinates'] * UnitLength * a * h**-1
    particles.velocity        = ggparticles['Velocity'] * UnitVelocity * a**(0.5)
    particles.mass            = ggparticles['Mass'] * UnitMass * h**-1
    particles.gadgetpotential = (ggparticles['Potential'] | units.cm * units.s**-1) * UnitVelocity * a**-1
    particles.gadgetID        = ggparticles['ParticleIDs']
    if ptype=="gas":
        particles.smoothinglength = ggparticles['SmoothingLength'] * UnitLength * a * h**-1
        particles.density         = ggparticles['Density'] * h**2 * a**-3 * UnitMass * UnitLength**-3
    return particles


if "__main__" in __name__:
    snip = h5py.File(snipfilename,'r')
    snipParameters = snip['Parameters']
    snipHeader = snip['Header']
    snipUnits = snip['Units']
    
    UnitMass = (snipUnits.attrs['UnitMass_in_g'] | units.g).as_quantity_in(units.MSun)
    UnitVelocity = (snipUnits.attrs['UnitVelocity_in_cm_per_s'] | units.cm * units.s**-1).as_quantity_in(units.kms)
    UnitTime = snipUnits.attrs['UnitTime_in_s'] | units.s
    UnitLength = (snipUnits.attrs['UnitLength_in_cm'] | units.cm).as_quantity_in(units.Mpc)
    
    h = snipHeader.attrs['HubbleParam']
    redshift = snipHeader.attrs['Redshift']
    a = (1./(1.+redshift))

    converter = nbody_system.nbody_to_si(UnitLength, UnitMass)
    
    #gasparticles  = readgadgethdfparticles(snip['PartType0'], ptype="gas")
    dmparticles   = readgadgethdfparticles(snip['PartType1'], ptype="dm")
    #dmparticles_a = readgadgethdfparticles(snip['PartType2'], ptype="dm")[0:15000]
    #dmparticles_b = readgadgethdfparticles(snip['PartType3'], ptype="dm")[0:15000]
    #starparticles = readgadgethdfparticles(snip['PartType4'], ptype="stars")

    print len(dmparticles)
    maxdist = 3 | units.Mpc
    loc = dmparticles[0].position#[60,157,197]|units.Mpc
    #print dmparticles[0].position
    maxdist2=maxdist**2

    dmparticles = dmparticles.select_array( lambda x,y,z: (x-loc[0])**2+(y-loc[1])**2+(z-loc[2])**2 < maxdist2,["x","y","z"] )
     
#    dmparticles = dmparticles.select(lambda x: abs(x-loc[0]) < maxdist, ['x'])
#    print len(dmparticles)
#    dmparticles = dmparticles.select(lambda x: abs(x-loc[1]) < maxdist, ['y'])
#    print len(dmparticles)
#    dmparticles = dmparticles.select(lambda x: abs(x-loc[2]) < maxdist, ['z'])
    print len(dmparticles)
#    dmparticles.distance2 = (dmparticles.x - loc[0])**2 + (dmparticles.y - loc[1])**2 + (dmparticles.z - loc[2])**2
#    dmparticles = dmparticles.select(lambda r2: r2 < maxdist**2, ['distance2'])
#    print len(dmparticles)

    raise Exception("hier")
    gravity = ph4(converter)

    #gravity.particles.add_particles(gasparticles)
    gravity.particles.add_particles(dmparticles)
    #gravity.particles.add_particles(dmparticles_a)
    #gravity.particles.add_particles(dmparticles_b)
    #gravity.particles.add_particles(starparticles)
    gravity.commit_particles()

    pos0 = dmparticles[0].position
#    pos1 = dmparticles[7001].position
#    pos2 = dmparticles[7002].position
#    pos0 = starparticles[0].position
#    pos1 = starparticles[1].position
#    pos2 = starparticles[2].position

    x, y, z = pos0
    pot_pos0 = gravity.get_potential_at_point(0.1|units.parsec,x,y,z)
#    x, y, z = pos1
#    pot_pos1 = gravity.get_potential_at_point(0.1|units.parsec,x,y,z)
#    x, y, z = pos2
#    pot_pos2 = gravity.get_potential_at_point(0.1|units.parsec,x,y,z)
    #pot_pos1 = gravity.get_potential_at_point(pos1)
    #pot_pos2 = gravity.get_potential_at_point(pos2)

    print dmparticles[0].gadgetpotential, pot_pos0, pot_pos0/dmparticles[0].gadgetpotential
#    print dmparticles_a[7000].gadgetpotential, pot_pos0, pot_pos0/dmparticles_a[7000].gadgetpotential
#    print dmparticles_a[7001].gadgetpotential, pot_pos1, pot_pos1/dmparticles_a[7001].gadgetpotential
#    print dmparticles_a[7002].gadgetpotential, pot_pos2, pot_pos2/dmparticles_a[7002].gadgetpotential
    #print starparticles[0].gadgetpotential, pot_pos0
    #print starparticles[1].gadgetpotential, pot_pos1
    #print starparticles[2].gadgetpotential, pot_pos2
