from job_server import JobServer

def somework(x):
  return x*x

def example_parallel_jobs(N):

  from socket import gethostname
  
  jobserver=JobServer([gethostname()])
  
  for i in range(N):
    jobserver.submit_job(somework, (i,))
   
  while jobserver.wait():
    job=jobserver.last_finished_job
    print job.args[0],job.result
    
if __name__=="__main__":
# this is needed in order for all the functions and classes to be 
# pickled with full module name
  from job_server_example import example_parallel_jobs  
  example_parallel_jobs(10)
