import os
import numpy

import matplotlib
from matplotlib import pyplot,gridspec,rc

from amuse.units import units
from amuse.datamodel import Particles
from amuse.io import read_set_from_file

c1 = '#e41a1c'
c2 = '#377eb8'
c3 = '#fc9272'
c4 = '#cccccc'

pi_180 = numpy.pi/180.

#from matplotlib.colors import LogNorm
rc('font',**{'size': 10})
#pyplot.rcParams['mathtext.fontset'] = "stix"
import brewer2mpl

def plot_all_snaps(bi, snapshot_dir, file_out,
                   center=01,
                   xmin=-250, xmax=250,
                   ymin=-250, ymax=250,
                   zmin=-50, zmax=120,
                   a1=-60.):
  
  try:
     os.mkdir(snapshot_dir)
  except:
     print "\t(" , snapshot_dir, "already exists)"
  
  fig = pyplot.figure(figsize=(6,8))
  lm = 0.1
  rm = 0.99
  tm = 0.95
  bm = 0.1
  gs1 = gridspec.GridSpec(2, 1, height_ratios=[ymax-ymin, zmax-zmin])
  gs1.update(left=lm, right=rm, top=tm, bottom=bm)
  xy_plane = pyplot.subplot(gs1[0, 0], aspect='equal')
  xz_plane = pyplot.subplot(gs1[1, 0], aspect='equal')
  
  # set center of the plot
  if center=='0':
    bi.position -= bi[0].position
    bi.velocity -= bi[0].velocity
  elif center=='1':
    bi.position -= bi[1].position
    bi.velocity -= bi[1].velocity
  elif center=='01':
    bi.move_to_center()
  
  n = len(bi)
  n_d = (n-2)/2
  i_d = (n-2)/2 + 2
  #filter_a = (bi.id == 'a')
  #filter_b = (bi.id == 'b')
  
  xb = bi[0:2].x.value_in(units.kpc)
  yb = bi[0:2].y.value_in(units.kpc)
  zb = bi[0:2].z.value_in(units.kpc)
  
  xd = bi[2:].x.value_in(units.kpc)
  yd = bi[2:].y.value_in(units.kpc)
  zd = bi[2:].z.value_in(units.kpc)
  
  xd_a = bi[2:i_d].x.value_in(units.kpc)
  yd_a = bi[2:i_d].y.value_in(units.kpc)
  zd_a = bi[2:i_d].z.value_in(units.kpc)
  xd_b = bi[i_d:].x.value_in(units.kpc)
  yd_b = bi[i_d:].y.value_in(units.kpc)
  zd_b = bi[i_d:].z.value_in(units.kpc)
  
  # rotate aroud z-axis
  a = a1*pi_180
  rot_z = ([numpy.cos(a), -numpy.sin(a), 0.0], 
            [numpy.sin(a), numpy.cos(a), 0.0], 
            [0.0, 0.0, 1.0])
  rot_a = (numpy.dot(rot_z, [xd_a.T, yd_a.T, zd_a.T])).T
  xd_a = rot_a[:,0]; yd_a = rot_a[:,1]; zd_a = rot_a[:,2]
  rot_a = (numpy.dot(rot_z, [xd_b.T, yd_b.T, zd_b.T])).T
  xd_b = rot_a[:,0]; yd_b = rot_a[:,1]; zd_b = rot_a[:,2]
  rot_a = (numpy.dot(rot_z, [xb.T, yb.T, zb.T])).T
  xb = rot_a[:,0]; yb = rot_a[:,1]; zb = rot_a[:,2]
  
  # plot name
  plot_i = snapshot_dir+file_out
  
  print " ** plotting: -> ", plot_i
  
  # XY plane
  xy_plane.scatter(xd_a, yd_a, c=c2, edgecolors = "none", alpha = 0.2, zorder=1)
  xy_plane.scatter(xd_b, yd_b, c=c3, edgecolors = "none", alpha = 0.2, zorder=1)
  xy_plane.scatter(xb[0], yb[0], c=c2, lw=0.5, zorder=2)
  xy_plane.scatter(xb[1], yb[1], c=c3, lw=0.5, zorder=2)
  
  # XZ plane
  xz_plane.scatter(xd_a, zd_a, c=c2, edgecolors = "none", alpha = 0.2, zorder=1)
  xz_plane.scatter(xd_b, zd_b, c=c3, edgecolors = "none", alpha = 0.2, zorder=1)
  xz_plane.scatter(xb[0], zb[0], c=c2, lw=0.5, zorder=2)
  xz_plane.scatter(xb[1], zb[1], c=c3, lw=0.5, zorder=2)
  
  # limits
  xy_plane.set_xlim(xmin, xmax)
  xy_plane.set_ylim(ymin, ymax)
  xz_plane.set_xlim(xmin, xmax)
  xz_plane.set_ylim(zmin, zmax)
  
  # labels
  xy_plane.set_xlabel('x [kpc]')
  xy_plane.set_ylabel('y [kpc]')
  xz_plane.set_xlabel('x [kpc]')
  xz_plane.set_ylabel('z [kpc]')
  
  #pyplot.tight_layout()
  #pyplot.savefig(plot_i, dpi=350)
  pyplot.savefig(plot_i)
  pyplot.cla()
    
  return
  