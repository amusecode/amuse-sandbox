#!/usr/bin/env python

import numpy, random, sys, time
from cap_w_mc_02 import call_cap_w
import emcee
from optparse import OptionParser
import socket
from amuse.ext.job_server import JobServer

#### kwistbeek setup
#sys.path.insert(0, '/data2/jilkova/amuse/kbo_flyby/')
#### para setup
#sys.path.insert(0, '/net/para35/data2/jilkova/emcee/dfm-emcee-c2e059f/')
#### LGM setup
#sys.path.insert(0, '/scratch/jilkova/amuse/kbo_flyby/') 
#sys.path.append('/scratch/jilkova/software/scipy/lib/python2.7/site-packages/') # scipy installation
#from emcee.utils import MPIPool
#from scipy import stats

begin_time = time.time()

### pool provided by the amuse job server
class JobServerPool(JobServer):
  def map(self,function,tasks):
    
    print 'tasks', len(tasks)
    js_out = open('job_server_map.out', 'a')
    print >> js_out, len(tasks), (time.time()-begin_time)/60.0/60.0
    js_out.close()
    
    result=[None]*len(tasks)
    for i,t in enumerate(tasks):
      job=self.submit_job(function, (t,))
      job.rank=i
    self.waitall()
    for job in self.finished_jobs:
      result[job.rank]=job.result
    return result

### sorting function
# sort the initial conditions for each iteration so
# that the fast walkers are calculated last
# !!!! constrains need to be updated by hand here !!!!
def sort_on_runtime(pos):
    p = numpy.atleast_2d(pos)
    idx_sort = numpy.empty((len(p[:,0])))
    constraints = numpy.array([0.2, 2.0,
                               200., 400.,
                               1., 4.0,
                               0., 360.,
                               0., 180.,
                               50., 200.])
    for i in xrange(len(p[:,0])):
      if not ( (constraints[0] < p[i,0] and p[i,0] < constraints[1]) and
               (constraints[2] < p[i,1] and p[i,1] < constraints[3]) and
               (constraints[4] < p[i,2] and p[i,2] < constraints[5]) and
               (constraints[6] < p[i,3] and p[i,3] < constraints[7]) and
               (constraints[8] < p[i,4] and p[i,4] < constraints[9]) and
              (constraints[10] < p[i,5] and p[i,5] < constraints[11]) ):
        # parameters outside the range => no flyby calculation, to do last
        #idx_sort[i] = numpy.inf
        idx_sort[i] = 1.
      else:
        # rest stays the same, i.e. random
        #idx_sort[i] = ((1.+p[i,2])*(1.+p[i,0]))/p[i,1]
        idx_sort[i] = 0.
    
    # fast jobs last
    idx = numpy.argsort(idx_sort[:])
    #print idx_sort
    #print p[idx]
    return p[idx], idx

### flyby calculation
def flyby(M = 1.0, Rp = 10.0, e = 1.0, i = 90.0, om = 90.0, rout = 100.0):
    
    rank, dist, sum_w_sed, sum_w_ased, nb, n_tot, rel_pos, rel_vel = \
      call_cap_w(m0=M,
                 peri=Rp,
                 ecc=e,
                 incl=i,
                 omega=om,
                 r_out=rout)
    
    print ' M', M, 'Rp', Rp, 'e', e, 'i', i, 'omega', om, 'rout', rout, 'rank', rank, 'dist', dist, 'w_sed', sum_w_sed, 'w_ased', sum_w_ased, 'nb', nb, 'ntot', n_tot, socket.gethostname()
    
    return rank, dist, rel_pos, rel_vel, sum_w_sed, sum_w_ased, nb, n_tot

# PROBABILITY FUNCTION used in EMCEE
def lnprobfn(x, constraints, f_output, directory, fname_i):
    
    ### check for constrains
    if not ( (constraints[0] < x[0] and x[0] < constraints[1]) and
             (constraints[2] < x[1] and x[1] < constraints[3]) and
             (constraints[4] < x[2] and x[2] < constraints[5]) and
             (constraints[6] < x[3] and x[3] < constraints[7]) and
             (constraints[8] < x[4] and x[4] < constraints[9]) and
            (constraints[10] < x[5] and x[5] < constraints[11]) ):
     
     rank = 0.
     dist = numpy.inf
     sum_w_sed = 0.
     sum_w_ased = 0.
     nb = 0
     n_tot = 0
     lnprob_final = -numpy.inf
     i_mco = -1
        
    ### if all parameters within constrains
    else:
      
      ### Running kuiperbelt
      try:
        rank, dist, rel_pos, rel_vel, sum_w_sed, sum_w_ased, nb, n_tot = \
          flyby(x[0], x[1], x[2], x[3], x[4], x[5])
        
        # in case no sednitos are transferred
        if (1. <= dist):
          lnprob_final = -numpy.inf
          i_mco = -2
        # in case some sednitos are tranferred
        else:
          # ln(1.-distance) 
          #   distance = 0 <=> 2d histograms are exactly the same
          #   distance = 1 <=> no overlap in 2d histograms
          lnprob_final = numpy.log(1.-dist)
          
          # to save the parameters and final relative positions and velocity vectors
          f_i = directory+fname_i+".dat"
          i_mco = numpy.loadtxt(f_i)
          i_mc = numpy.array([i_mco+1])
          numpy.savetxt(f_i, i_mc)
          f_out_i = directory+"run_{0:05d}.npy".format(int(i_mco))
          f_out_i_par = directory+"run_{0:05d}.par".format(int(i_mco))
          numpy.save(f_out_i, [rel_pos, rel_vel])
          numpy.savetxt(f_out_i_par, x)
      
      ### in case of crash / error of flyby calculation -- should not happen often
      except:
        print " ** !!! exception -- error in flyby !!!"
        lnprob_final = -numpy.inf
        rank = 0.
        dist = numpy.inf
        sum_w_sed = 0.
        sum_w_ased = 0.
        nb = 0
        n_tot = 0
        i_mco = -666
    
    ### output to save
    out_p = [x[0], x[1], x[2], x[3], x[4], x[5], n_tot, nb, sum_w_sed, sum_w_ased, rank, dist, lnprob_final, int(i_mco)]
    
    output = open(f_output, 'a')
    print >> output,  "{} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(*out_p)
    output.close()
    
    return lnprob_final

def main(file_name, directory, nwalkers, nburn, nchain, fout, nthreads,
         mmin, mmax, pmin, pmax, emin, emax, imin, imax, omin, omax, rmin, rmax,
         fname_i, f_input, f_ini):
    ndim = 6
    constraints = [mmin, mmax,   # M_min, M_max [MSun]
                   pmin, pmax,   # R_p_min, R_p_max [AU]
                   emin, emax,   # e_min, e_max
                   imin, imax,   # i_min, i_max [deg]
                   omin, omax,   # omega_min, omega_max [deg]
                   rmin, rmax]   # outer radius [AU]
    print ' ** adopted constrains:'
    print '\t', constraints[:]
    
    # to delate the input and output files at the beginning
    f_output = fout
    output = open(f_output, 'w')
    output.close()
    
    # create the file to count the runs and save 0
    f_i = directory+fname_i+".dat"
    i0 = numpy.array([0])
    numpy.savetxt(f_i, i0)
    
    # walkers initialization
    print ' ** Running MCMC to determine the initial conditions of the perturber + disk'
    print ' ** Initializing ensemble with ', str(ndim), ' dimensions and ', str(nwalkers), ' walkers'
    x_i = numpy.empty((nwalkers,ndim))

    # use some sucesful initial conditions from previous MC runs
    #   read from file f_ini if it exists
    try:
      x_ini = numpy.loadtxt(f_ini)
      n_old = len(x_ini[:,0])
      x_i[-n_old:] = x_ini
      print " \t\t Reading walkers from", f_ini, ",", n_old, "walkers"
    except:
      n_old = 0
      
    M_i  = [numpy.random.uniform(constraints[0], constraints[1]) for j in xrange(nwalkers-n_old)]
    Rp_i = [numpy.random.uniform(constraints[2], constraints[3]) for j in xrange(nwalkers-n_old)]
    e_i  = [numpy.random.uniform(constraints[4], constraints[5]) for j in xrange(nwalkers-n_old)]
    i_i  = [numpy.random.uniform(constraints[6], constraints[7]) for j in xrange(nwalkers-n_old)]
    om_i = [numpy.random.uniform(constraints[8], constraints[9]) for j in xrange(nwalkers-n_old)]
    ro_i = [numpy.random.uniform(constraints[10], constraints[11]) for j in xrange(nwalkers-n_old)]
    
    x_i[:nwalkers-n_old] = [numpy.array([M_i[j], Rp_i[j], e_i[j], i_i[j], om_i[j], ro_i[j]],) \
      for j in xrange(nwalkers-n_old)]
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprobfn, 
                                    args=[constraints, f_output, directory, fname_i],
                                    pool=pool, 
                                    runtime_sortingfn=sort_on_runtime)
    
    print ' ** Entering "burn in" phase with ', str(nburn), ' burn in iterations'
    pos, prob, state = sampler.run_mcmc(x_i, nburn)
    input = open(f_input, 'a')
    print >> input, pos, prob, state
    input.close()
    sampler.reset()
    
    print ' ** Running detailed MCMC from burned in positions with ', str(nchain), ' iterations..\n'
    pos2, prob2, state2 = sampler.run_mcmc(pos, nchain, rstate0=state)
    input = open(f_input, 'a')
    print >> input, pos2, prob2, state2
    input.close()
    
    print "Mean acceptance fraction:", numpy.mean(sampler.acceptance_fraction)
    print "Acceptance fractions:", sampler.acceptance_fraction
    print '\nInput written out in ', f_input
    print 'Output written out in ', f_output
    
    run_time = time.time() - begin_time
    print '\n', "Program finished successfully!"
    if run_time/60 < 1.0: print "Run time: " + str(run_time) + " sec.", '\n'
    else:                 print "Run time: " + str(run_time/60) + " min.", '\n'

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="file_name", type="str", default='sednas_stones.npy', 
                      help="The observed KBO to reproduce")
    result.add_option("--dir", dest="directory", type="str", default='./', 
                      help="The working directory")
    result.add_option("--nw", dest="nwalkers", type="int", default=100, 
                      help="Number of walkers of MCMC run")
    result.add_option("--nb", dest="nburn", type="int", default=100, 
                      help="Number of burn in iterations of MCMC run")
    result.add_option("--nc", dest="nchain", type="int", default=1000,
                      help="Number of detailed iteratoins of MCMC run")
    result.add_option("--fout", dest="fout", type="str", default='mcmc_out.dat', 
                      help="file name to write results")
    result.add_option("--f_input", dest="f_input", type="str", default='mcmc_sampler.dat', 
                      help="file name to write sampler after burn-in and chain phases")
    result.add_option("--f_ini", dest="f_ini", type="str", default='mcmc_ini.dat', 
                      help="file name with initial conditions for burn-in")
    result.add_option("--nt", dest="nthreads", type="int", default=1, 
                      help="number of threads for emcee")
    result.add_option("--mmin", dest="mmin", type="float", default = 0.02,
                      help="perturbers minimal mass in MSun [0.3]")
    result.add_option("--mmax", dest="mmax", type="float", default = 0.05,
                      help="perturbers maximal mass in MSun [3.0]")
    result.add_option("--pmin", dest="pmin", type="float", default = 40.0,
                      help="perturbers minimal pericenter in AU [40]")
    result.add_option("--pmax", dest="pmax", type="float", default = 300.0,
                      help="perturbers maximal pericenter in AU [1000]")
    result.add_option("--emin", dest="emin", type="float", default = 1.0,
                      help="perturberts minimal eccentricity [1.0]")
    result.add_option("--emax", dest="emax", type="float", default = 5.0,
                      help="perturbers maximal eccentricity [100.0]")
    result.add_option("--imin", dest="imin", type="float", default = 0.0,
                      help="perturbers minimal inclination in deg [0.0]")
    result.add_option("--imax", dest="imax", type="float", default = 180.0,
                      help="perturbers maximal inclination in deg [180.0]")
    result.add_option("--omin", dest="omin", type="float", default = 0.0,
                      help="perturbers minimal argument of pericenter in deg [0.0]")
    result.add_option("--omax", dest="omax", type="float", default = 180.0,
                      help="perturbers maximal argument of pericenter in deg [180.0]")
    result.add_option("--rmin", dest="rmin", type="float", default = 50.0,
                      help="perturbers minimal argument of pericenter in deg [0.0]")
    result.add_option("--rmax", dest="rmax", type="float", default = 200.0,
                      help="perturbers maximal argument of pericenter in deg [180.0]")
    result.add_option("--fname_i", dest="fname_i", type="string", default = "i_mc",
                      help="filename to store the runs ID")
    return result

if __name__ in ('__main__'):
  
  ### set up distributed amuse
  import mc_distributed
  
  # set up nodes (only names)
  #lgm_node_names = ["node02", "node03"]
  lgm_node_names = ["node17", "node15", "node13", "node12", "node11", "node09", "node08", "node07", "node06", "node04", "node03"]
  
  # initialize distributed amuse
  # number of CPUs and times have to be specified in mc_distributed
  instance=mc_distributed.start_distributed(lgm_node_names=lgm_node_names)
  
  ### set up job server
  # set up all host CPUs = nodes * number_of_cpus
  #hosts=["node02"]*10+["node03"]*10 
  hosts=["node17"]*12+["node15"]*12+["node13"]*12+["node12"]*12+["node11"]*12+ \
        ["node09"]*12+["node08"]*12+["node07"]*12+["node06"]*12+["node04"]*12+ \
        ["node16"]*12+["node03"]*12

  # no_wait = False ---> workers __do__ wait
  from mcmc_kuiperbelt_10 import lnprobfn
  pool=JobServerPool(hosts, channel_type="distributed", no_wait=False)
  
  ### get arguments
  o, arguments  = new_option_parser().parse_args()
  
  ### run emcee, run!
  main(**o.__dict__)
  
  ### time for you to leave ###
