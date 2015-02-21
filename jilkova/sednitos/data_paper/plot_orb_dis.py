"""
Data and plots for Sednitos paper
"""

import numpy
import argparse

from matplotlib import pyplot, rc, gridspec
from matplotlib.colors import LogNorm

rc('text', usetex=True)
rc('font',**{'family':'sans-serif', 'size': 12})

pi_180 = numpy.pi/180.0

### colors
blue = '#4393c3'
red = '#f03b20'
almost_black = '#262626'

def multiplot_var_par2(parx, par1, par2,
                       plot_name='plot_orb_dis',
                       lx='', ly1='', ly2='',
                       minx=None, maxx=None, 
                       miny1=None, maxy1=None,
                       miny2=None, maxy2=None,
                       star_ini=None, 
                       log=False,
                       plot_each=1,
                       tit=None):
  """
  Fig.2
  -- par1 and par2 vs. parx
  """
  
  ### markers
  msize=6     # size of the points
  c_e='none'  # edge_color
  l_o=0.0     # edge width
  a=0.5       # trnasparency
  
  ###
  #msize=30
  #a=0.2
  #l_o=0.4
  #c_e=almost_black
  
  ### set multiplot
  fig = pyplot.figure(figsize=(5,8))
  gs1 = gridspec.GridSpec(2, 1)
  ax = pyplot.gca()
  lm = 0.13
  rm = 0.97
  tm = 0.96
  bm = 0.06
  hs = 0.0
  ws = 0.0
  gs1.update(left=lm, right=rm, top=tm, bottom=bm, wspace=ws, hspace=hs)
  pl2 = pyplot.subplot(gs1[0,0])
  pl1 = pyplot.subplot(gs1[1,0])
  
  ### set log-scale
  if log is True:
    pl1.set_xscale('log')
    pl2.set_xscale('log')
  
  ### filters
  filter0 = ( star_ini==0 )  # initially Sun's stone
  filter1 = ( star_ini==1 )  # initially Q's stone
  
  ### set order and colors for plotting
  # fixed -- Sun red and on top
  filters = [filter1, filter0]
  colors = [blue, red]
  
  # adjusted -- star with more particles on top
  #if sum(filter0) > sum(filter1):
    #filters = [filter0, filter1]
    #colors = [red, blue]
  #else:
    #filters = [filter1, filter0]
    #colors = [blue, red]
  
  ### plotting
  st1=pl1.scatter(parx[filters[0]][::plot_each], par1[filters[0]][::plot_each], 
              alpha=a, marker='o', 
              edgecolors=c_e, 
              c=colors[0], 
              s=msize,
              linewidth=l_o,
              rasterized=True)
  st2=pl1.scatter(parx[filters[1]][::plot_each], par1[filters[1]][::plot_each], 
              alpha=a, marker='o', 
              edgecolors=c_e, 
              c=colors[1], 
              s=msize,
              linewidth=l_o,
              rasterized=True)
  
  pl2.scatter(parx[filters[0]][::plot_each], par2[filters[0]][::plot_each], 
              alpha=a, marker='o', 
              edgecolors=c_e, 
              c=colors[0], 
              s=msize,
              linewidth=l_o,
              rasterized=True)
  pl2.scatter(parx[filters[1]][::plot_each], par2[filters[1]][::plot_each], 
              alpha=a, marker='o', 
              edgecolors=c_e, 
              c=colors[1], 
              s=msize,
              linewidth=l_o,
              rasterized=True)
      
  ### title
  pyplot.text(0.5, 1.08, tit, 
              horizontalalignment='center', verticalalignment='bottom', 
              transform=ax.transAxes)
    
    
  ### counts of bound particles
  #bound_ini0 = sum(filter0)
  #bound_ini1 = sum(filter1)
  #pyplot.text(0.95, 1.08, r'${0:s}$'.format(str(bound_ini0))+' /',
              #horizontalalignment='right', verticalalignment='bottom', 
              #transform=ax.transAxes)
  #pyplot.text(1.06, 1.08, r'${0:s}$'.format(str(bound_ini1))    ,
              #horizontalalignment='right', verticalalignment='bottom', 
              #transform=ax.transAxes)
    
  ### setting labels and limits -- should be done in more pythonic way
  pl1.set_xlabel(lx)
  pl2.set_xlabel('')
  pl1.set_ylabel(ly1)
  pl2.set_ylabel(ly2)
  pyplot.setp(pl2.get_xticklabels(), visible=False)
  #pl1.set_xlim(xl_min1, xl_max1)
  #pl2.set_xlim(xl_min2, xl_max2)
  #pl1.set_ylim(yl_min1, yl_max1)
  #pl2.set_ylim(yl_min2, yl_max2)
  pl1.set_xlim(minx, maxx)
  pl2.set_xlim(minx, maxx)
  pl1.set_ylim(miny1, maxy1)
  pl2.set_ylim(miny2, maxy2)
  
  ### save plot as .png
  plot_png = plot_name+".png"
  print " ** plot", " -> ", plot_png
  pyplot.savefig(plot_png)
  
  ### save plot as .pdf
  #plot_pdf = plot_name+".pdf"
  #pyplot.savefig(plot_pdf, dpi=400)
  
  return
  
def plot_2orpar_rini(par1x, par1y, par2y, par_color,
                     par1x_o, par1y_o, par2y_o,
                     lx='', ly1='', ly2='', lc='',
                     plot_out='plot_orb_par_rini',
                     minx=None, maxx=None, 
                     miny1=None, maxy1=None,
                     miny2=None, maxy2=None,
                     plot_each=1):
  """
  Fig.1
  -- orbital distribution color mapped by initial radius
  -- par1 and par2 vs. parx
  -- simulated and observed data
  """

  ### set multiplot
  fig = pyplot.figure(figsize=(6,8))
  gs1 = gridspec.GridSpec(2, 1)
  ax = pyplot.gca()
  lm = 0.1
  rm = 0.83
  tm = 0.98
  bm = 0.06
  hs = 0.0
  ws = 0.0
  gs1.update(left=lm, right=rm, top=tm, bottom=bm, wspace=ws, hspace=hs)
  pl2 = pyplot.subplot(gs1[0,0])
  pl1 = pyplot.subplot(gs1[1,0])
  
  ### color map
  #cmap = pyplot.get_cmap('YlOrRd', 6)
  cmap = pyplot.get_cmap('YlGnBu',10)
  #cmap = cm.get_cmap('cubehelix_r', 8)
  
  ### markers
  c_o = red
  #c_o = 'none'
  m_o = 'D'
  s_o = 70
  l_o = 2
  l_s=0.1
  #e_s = 'gainsboro'
  e_s = almost_black
  e_o = almost_black
  a_o=0.8
  
  ### plots
  s1=pl1.scatter(par1x, par1y, c=par_color, cmap=cmap,
                 edgecolors=e_s,
                 linewidth=l_s, 
                 s=30,
                 rasterized=True)
  pl1.scatter(par1x_o, par1y_o, 
              facecolor=c_o,
              marker=m_o,
              s=s_o,
              edgecolors=e_o,
              linewidth=l_o,
              alpha=a_o)
  
  s2=pl2.scatter(par1x, par2y, c=par_color, cmap=cmap,
                 edgecolors=e_s,
                 linewidth=l_s, 
                 s=30,
                 rasterized=True)
  pl2.scatter(par1x_o, par2y_o, 
              facecolor=c_o,
              marker=m_o,
              s=s_o,
              edgecolors=e_o,
              linewidth=l_o,
              alpha=a_o)
  ###

  ### set labels and limits
  pl1.set_xlabel(lx)
  pl2.set_xlabel('')
  pl1.set_ylabel(ly1)
  pl2.set_ylabel(ly2)
  pl1.set_xlim(minx, maxx)
  pl2.set_xlim(minx, maxx)
  pl1.set_ylim(miny1, maxy1)
  pl2.set_ylim(miny2, maxy2)
  pyplot.setp(pl2.get_xticklabels(), visible=False)
  
  ### colorbar
  cbaxes = fig.add_axes([rm+0.03, bm, 0.03, (tm-bm)/2. ])
  cbar = pyplot.colorbar(s1, cax=cbaxes)
  cbar.set_label(lc)
  
    ### save plot as .png
  plot_png = plot_out+".png"
  print " ** plot -> ", plot_png
  pyplot.savefig(plot_out)
  
  ### save plot as .pdf
  #plot_pdf=plot_out+'.pdf'
  #print " ** plot -> ", plot_pdf
  #pyplot.savefig(plot_pdf, dpi=300)
  
  return
  
def transform_m180_p180(angles):
  angles_tr = angles
  for i,a_i in enumerate(angles):
    if (a_i>180.0):
      angles_tr[i] = a_i-360.0
  return angles_tr

def new_option_parser():
  result = argparse.ArgumentParser()
  result.add_argument("--fin", 
                    dest="file_in",  default="run_t100_disk.npz", type=str,
                    help="input file")
  result.add_argument("--fout", 
                    dest="file_out",  default="plot_orb_dis", type=str,
                    help="output file")
  result.add_argument("--fobs", 
                    dest="file_obs",  default="sednas_stones.npy", type=str,
                    help="observed sednitos")
  result.add_argument("--amin", 
                    dest="a_min", default = 0,
                    help="minimal semi-major axis")
  result.add_argument("--amax", 
                    dest="a_max", default = 700,
                    help="maximal semi-major axis")
  result.add_argument("--emin", 
                    dest="e_min", default = 0,
                    help="minimal eccentricity")
  result.add_argument("--emax", 
                    dest="e_max", default = 1,
                    help="maximal eccentricity")
  result.add_argument("--qmin", 
                    dest="q_min", default = None,
                    help="minimal perihelion")
  result.add_argument("--qmax", 
                    dest="q_max", default = None,
                    help="maximal perihelion")
  result.add_argument("--imin", 
                    dest="i_min", default = -10,
                    help="minimal inclination")
  result.add_argument("--imax", 
                    dest="i_max", default = 100,
                    help="maximal inclination")
  result.add_argument("--f_star_ini", 
                    dest="f_star_ini", default="star_ini.npy", type=str,
                    help="file with initial disk membership")
  result.add_argument("--f_sem_ini", 
                    dest="f_sem_ini", default="sem_ini.npy", type=str,
                    help="file with initial disk radius")
  result.add_argument("--r0_max", 
                    dest="r0_max", default=90.0, type=float,
                    help="max radius for the Sun's disk")
  result.add_argument("--plot_each", 
                    dest="plot_each", default=1, type=int,
                    help="plot each plot_each point")
  return result

if __name__ in ('__main__', '__plot__'):
  arg = new_option_parser().parse_args()
  
  ### get the initial membership (0 -- Sun, 1 -- Q)
  star_ini = numpy.load(arg.f_star_ini)
  
  ### get the initial semi-major axis
  sem_ini = numpy.load(arg.f_sem_ini)
  
  ### get the final orbital elements from the .npz file
  # ..._t0 -- with respect to the Sun
  # ..._t1 -- with respect to Q
  orb_par = numpy.load(arg.file_in)
  sem_t0=orb_par['sem_t0']          # semi-major axis, Sun, [AU]
  ecc_t0=orb_par['ecc_t0']          # eccentricity, Sun
  per_t0=orb_par['per_t0']          # period, Sun, [yr]
  inc_t0=orb_par['inc_t0']/pi_180   # inclination, Sun, [deg]
  arg_t0=orb_par['arg_t0']/pi_180   # argument of perihelion, Sun, [deg]
  lon_t0=orb_par['lon_t0']/pi_180   # longitude of ascending node, Sun, [deg]
  sem_t1=orb_par['sem_t1']          # and the smae for Q ...
  ecc_t1=orb_par['ecc_t1']
  per_t1=orb_par['per_t1']
  inc_t1=orb_par['inc_t1']/pi_180
  arg_t1=orb_par['arg_t1']/pi_180
  lon_t1=orb_par['lon_t1']/pi_180
  orb_par.close()
  
  ### filter by initial radius if r0_max for the Sun's disk is defined
  if (arg.f_sem_ini is not None) and (arg.r0_max is not None):
    filter_r0 = ( ((star_ini==0) & (sem_ini<=arg.r0_max)) | (star_ini==1) )
    print " ** filtering disk around the Sun, R_max =", arg.r0_max, "AU"
    print " \t # of stones in Sun's disk", sum(filter_r0)-sum(star_ini)
    sem_t0 = sem_t0[filter_r0]
    sem_t1 = sem_t1[filter_r0]
    ecc_t0 = ecc_t0[filter_r0]
    ecc_t1 = ecc_t1[filter_r0]
    per_t0 = per_t0[filter_r0]
    per_t1 = per_t1[filter_r0]
    inc_t0 = inc_t0[filter_r0]
    inc_t1 = inc_t1[filter_r0]
    arg_t0 = arg_t0[filter_r0]
    arg_t1 = arg_t1[filter_r0]
    lon_t0 = lon_t0[filter_r0]
    lon_t1 = lon_t1[filter_r0]
    star_ini = star_ini[filter_r0]
    sem_ini = sem_ini[filter_r0]
    
  ### get perihelion distance
  q_t1 = (numpy.ones_like(ecc_t1)-ecc_t1)*sem_t1
  q_t0 = (numpy.ones_like(ecc_t0)-ecc_t0)*sem_t0
  
  ### read observed data
  # sednitos
  orb_par_obs = numpy.load(arg.file_obs)
  sem_o = orb_par_obs[:,2]
  ecc_o = orb_par_obs[:,3]
  q_o = orb_par_obs[:,0]
  arg_o_360 = orb_par_obs[:,5]
  arg_o = transform_m180_p180(arg_o_360)
  inc_o_360 = orb_par_obs[:,4]
  inc_o = transform_m180_p180(inc_o_360)
  
  ### plotting 
  
  ### Fig. 1
  # filter: initially Q's disk, bound to the Sun
  filter_01 = ( (star_ini==1) & (ecc_t0[:,0]<1.0) )
  plot_out = arg.file_out+'_rini'
  plot_2orpar_rini(sem_t0[(filter_01),0], 
                   ecc_t0[(filter_01),0], 
                   inc_t0[(filter_01),0], 
                   sem_ini[filter_01],
                   sem_o, ecc_o, inc_o,
                   lx=r'$a\,[\mathrm{AU}]$', 
                   ly1=r'$e$',
                   ly2=r'$i\,[\mathrm{deg}]$',
                   lc=r'$r_{\mathrm{ini,Q}}\,[\mathrm{AU}]$',
                   plot_out=plot_out,
                   minx=arg.a_min, maxx=arg.a_max, 
                   #miny1=arg.e_min, 
                   maxy1=arg.e_max,
                   miny2=arg.i_min, maxy2=arg.i_max,
                   plot_each=1)
  
  ### Fig. 2
  # Q's disk
  # to filter stones bound to Q
  filter_1 = ( ecc_t1<1. )
  plot_name = arg.file_out+'_a1'
  multiplot_var_par2(sem_t1[filter_1], ecc_t1[filter_1], inc_t1[filter_1],
                     plot_name=plot_name,
                     lx=r'$a\,[\mathrm{AU}]$', 
                     ly1=r'$e$',
                     ly2=r'$i\,[\mathrm{deg}]$',
                     minx=arg.a_min, maxx=arg.a_max, 
                     miny1=arg.e_min, maxy1=arg.e_max,
                     miny2=arg.i_min, maxy2=arg.i_max,
                     star_ini=star_ini[filter_1[:,0]],
                     plot_each=arg.plot_each,
                     tit=r'$\mathrm{Q:\,star\,\,of\,\,the\,\,preferred\,\,encounter}$')
                     
  # Sun's disk
  # to filter stones bound to the Sun
  filter_0 = ( ecc_t0<1. )
  plot_name = arg.file_out+'_a0'
  multiplot_var_par2(sem_t0[filter_0], ecc_t0[filter_0], inc_t0[filter_0], 
                     plot_name=plot_name,
                     lx=r'$a\,[\mathrm{AU}]$', 
                     ly1=r'$e$',
                     ly2=r'$i\,[\mathrm{deg}]$',
                     minx=arg.a_min, maxx=arg.a_max, 
                     miny1=arg.e_min, maxy1=arg.e_max,
                     miny2=arg.i_min, maxy2=arg.i_max,
                     star_ini=star_ini[filter_0[:,0]],
                     plot_each=arg.plot_each,
                     tit=r'$\mathrm{Sun}$')
  
  ###
  