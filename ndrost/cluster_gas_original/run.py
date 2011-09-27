import numpy
import clustergas_gadget as clustergas
from amuse.support.units import units

numpy.random.seed(123477)

clustergas.clustergas_restart(
                   "cl_Ns1k_Ng100k_sf03_sn01_R05_s1",
                   238)

#clustergas.clustergas(sfeff=0.3, 
#                    Nstar=1000,
#                    Ngas=100000,
#                    Rscale=0.5 | units.parsec,
#                    feedback_efficiency=0.1,
#                    runid="cl_Ns1k_Ng100k_sf03_sn01_R05_s1")

#clustergas.merge_clusters_with_gas(Nstar=[5000,5000],totalNgas=400000,
#                    runid="mergetest2",tmerge=2. | units.Myr, 
#                    t_end=7.5 | units.Myr, dt_plot=0.025 | units.Myr,
#                    eps_star=0.00025 | units.parsec)

#clustergas.clustergas_restart(
#                   "cl_gdg_Ns1k_Ng100k_sf03_sn001_R01_s2",
#                   659)

#numpy.random.seed(123491)
#clustergas.clustergas(sfeff=0.05, 
#                   Nstar=1000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.01,
#                   runid="cl_gdg_Ns1k_Ng100k_sf005_sn001_R05_s2")


#numpy.random.seed(12351)
#clustergas.clustergas(sfeff=0.5, 
#                   Nstar=100,
#                   Ngas=10000,t_end=.35 | units.Myr,
#                   Rscale=0.25 | units.parsec,
#                   feedback_efficiency=0.1,
#                   runid="test")

#clustergas.clustergas(sfeff=0.5, 
#                   Nstar=1000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.1,
#                   runid="cl_gdg_Ns1k_Ng100k_sf05_sn01_R05_s2")

#clustergas.clustergas_restart(
#                   "cl_gdg_Ns3k_Ng100k_sf05_sn001_R05_s1",
#                   457)

#clustergas.clustergas(sfeff=0.5, 
#                   Nstar=3000,
#                   Ngas=300000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.01,
#                   runid="cl_gdg_Ns3k_Ng300k_sf05_sn001_R05_s1")

#clustergas.clustergas(sfeff=0.75, 
#                   Nstar=3000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.01,
#                   runid="cl_gdg_Ns1k_Ng100k_sf075_sn001_R05")

#clustergas.clustergas(sfeff=0.5, 
#                   Nstar=1000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.01,
#                   runid="cl_gdg_Ns1k_Ng100k_sf05_sn001_R05")

#clustergas.clustergas(sfeff=0.5, 
#                   Nstar=1000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.1,
#                   runid="cl_gdg_Ns1k_Ng100k_sf05_sn01_R05")

#clustergas.clustergas(sfeff=0.25, 
#                   Nstar=1000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.01,
#                   runid="cl_gdg_Ns1k_Ng100k_sf025_sn001_R05-2")

#clustergas.clustergas(sfeff=0.05, 
#                   Nstar=1000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.01,
#                   runid="cl_gdg_Ns1k_Ng100k_sf005_sn001_R05-2")


#clustergas.clustergas(sfeff=0.05, 
#                   Nstar=1000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.01,
#                   runid="cl_Ns1k_Ng100k_sf005_sn001_R05")

#clustergas.clustergas(sfeff=0.05, 
#                   Nstar=1000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.001,
#                   runid="cl_Ns1k_Ng100k_sf005_sn0001_R05")

#clustergas.clustergas(sfeff=0.25, 
#                   Nstar=1000,
#                   Ngas=100000,
#                   Rscale=0.5 | units.parsec,
#                   feedback_efficiency=0.001,
#                   runid="cl_gdg_Ns1k_Ng100k_sf025_sn0001_R05")


