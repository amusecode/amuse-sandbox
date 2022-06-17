import numpy

from matplotlib import pyplot


def regular_resample(Ntarget,ind):
    """ regular resampling """
    N=len(ind)
    if N<=Ntarget:
        return ind
    ind_=numpy.mgrid[0:N-1:Ntarget*1j]
    ind_=numpy.array(ind_,int)
    return ind[ind_]

def resample_tracking_extrema(Ntarget,ind,arr):
    """ resample while attempting to track extrema """
    if len(ind)<=Ntarget:
        return ind

    ind=regular_resample(Ntarget,ind)

    for ithis in range(1,len(ind)):
      ntarget=ind[ithis]
      ncurr=ind[ithis-1]
      while ncurr < ntarget-1:
        nprev=ncurr
        ncurr=ncurr+1
        nnext=ncurr+1
        if ((arr[ncurr]-arr[nprev])*(arr[nnext]-arr[ncurr])) < 0* arr[ncurr]:
          ind[ithis]=ncurr
    return ind
  
def resample_minimize(Ntarget,ind, x, *arr):
    """ resample minimizing the mean root square error of the linear interpolant of the input arrs 
    
    returns indexing of the resampling
    """ 
    N=len(ind)
    if N<=Ntarget:
        return ind

    err=0*arr[0][1:-1]
    for a in arr:
      intpol=a[:-2]+(x[1:-1]-x[:-2])*(a[2:]-a[:-2])/(x[2:]-x[:-2])
      err+=(a[1:-1]-intpol)**2

    index=numpy.argsort(err)+1
    
    Ndrop=N-Ntarget
    dropindex=index[:Ndrop]
    dropindex=sorted(dropindex)
    
    isolated=[]
    prev=-1
    for i in dropindex:
      if prev+1!=i:
        isolated.append(i)
        prev=i
    
    ind[isolated]=-1
    
    a=numpy.where(ind>=0)
    ind=ind[a]

    if len(ind)>Ntarget:
        ind_=resample_minimize(Ntarget, numpy.arange(len(ind)), x[ind], *[a[ind] for a in arr])
        ind=ind[ind_]
    
    return ind

if __name__=="__main__":
    N=100
    Ntarget=9
    ind=numpy.arange(N)
    x=numpy.arange(N)/N
    #~ y=(numpy.arange(N)/N)
    y=(numpy.arange(N)/N)**2
    #~ y=numpy.sin(10*numpy.arange(N)/N)
    res=regular_resample(Ntarget,ind)
    assert len(res)==Ntarget
    print(res)

    res=resample_tracking_extrema(Ntarget,ind,y)
    assert len(res)==Ntarget
    print(res)

    res=resample_minimize(Ntarget,ind,x, y)
    assert len(res)==Ntarget
    print(res)

