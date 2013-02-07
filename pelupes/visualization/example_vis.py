"""
 prototype for Maarten
"""

from amuse.visualization.Vis3D

from amuse.ic.plummer import new_plummer_model

plum=new_plummer_model(100)

vis=Vis3D()

plum.r=1
plum.g=1
plum.b=1

vis.particles.add_particles(plum)

vis.start_viewer()

vis.stop()
