from amuse.io import read_set_from_file,write_set_to_file
from amuse.units import units,generic_unit_converter
from amuse.datamodel import ParticlesWithUnitsConverted,Particles

s=read_set_from_file("example.dat","gadget")

conv=generic_unit_converter.ConvertBetweenGenericAndSiUnits(3.085678e21 | units.cm, 1.989e43 | units.g, 1e5 | units.cm / units.s)

p=ParticlesWithUnitsConverted(s[1],conv.as_converter_from_si_to_generic())
pp=Particles(len(p))
pp.mass=p.mass
pp.x=p.x
pp.y=p.y
pp.z=p.z
pp.vx=p.vx
pp.vy=p.vy
pp.vz=p.vz

write_set_to_file(pp, "example.amuse","amuse")

