from job_server import JobServer

def somework(x):
  return x*x

if __name__=="__main__":
  from socket import gethostname
  print gethostname()
  
  jobserver=JobServer([gethostname()]*4)
  
  for i in range(10):
    jobserver.submit_job(somework, (i,))
   
  while jobserver.wait():
    job=jobserver.last_finished_job
    print job.args[0],job.result
    
  
