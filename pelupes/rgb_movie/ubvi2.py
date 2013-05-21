import numpy

from amuse.units import units,constants,nbody_system

from amuse.io import read_set_from_file

from amuse.ext.job_server import JobServer

import Image

def calculate_effective_temperature(luminosity, radius):
    return ((luminosity/(constants.four_pi_stefan_boltzmann*radius**2))**.25).in_(units.K)

if __name__=="__main__":
  import time
  from ubvivariant import rgb_frame
  directory="../"
  label="sse"
  L=20. | units.parsec
  image_size=[1920,1080]
  percentile=0.999

  i=0
  stars=read_set_from_file(directory+label+"-se-%6.6i"%i,"amuse",close_file=True)
  grav=read_set_from_file(directory+label+"-grav-%6.6i"%i,"amuse",close_file=True)
  channel=grav.new_channel_to(stars)
  channel.copy_attributes(["x","y","z"])

  stars.temperature=calculate_effective_temperature(stars.luminosity,stars.radius)

  vmax=rgb_frame(stars,dryrun=True,image_width=L,image_size=image_size,percentile=percentile)
  
  t0=time.time()
  nfirst=0
  nsnap=1000
  
  jobserver=JobServer(channel_type="mpi",hosts=["emphyrio"]*4 )

  for i in range(nfirst,nsnap+1):
    stars=read_set_from_file(directory+label+"-se-%6.6i"%i,"amuse",close_file=True)
    grav=read_set_from_file(directory+label+"-grav-%6.6i"%i,"amuse",close_file=True)
    channel=grav.new_channel_to(stars)
    channel.copy_attributes(["x","y","z"])
  
    stars.temperature=calculate_effective_temperature(stars.luminosity,stars.radius)
    
    job=jobserver.submit_job(rgb_frame,
      args=(stars,), 
      kwargs=dict(dryrun=False, vmax=vmax, image_width=L, image_size=image_size, percentile=percentile) )
    job.i=i
  
  i=nfirst
  while jobserver.wait():
    i+=1
    job=jobserver.last_finished_job
    vmax,image=job.result
    imrgb=Image.fromstring(image['mode'], image['size'], image['pixels'])
    imrgb.save("./frames/"+label+"-rgb-%6.6i.png"%job.i)
      
    t2=time.time()
    eta=(nsnap-i)*(t2-t0)/(i-nfirst+1)
    print i, "eta (minutes):",eta/60.
