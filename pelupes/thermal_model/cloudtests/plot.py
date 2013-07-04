import numpy

from matplotlib import pyplot

from amuse.io import read_set_from_file

from amuse.units import units

from StringIO import StringIO

from amuse.ext.radial_profile import radial_density

def mf_histo(parts):
  
    masses = parts.mass.value_in(units.MSun)
    
    bins = 10**numpy.linspace(-2, 2, 20)
    number_of_particles, bin_edges= numpy.histogram(masses, bins = bins)
    
    bin_sizes = bin_edges[1:] - bin_edges[:-1]
    
    y = number_of_particles / bin_sizes
#    x = (bin_edges[1:] + bin_edges[:-1]) / 2.0

    xx=[]
    yy=[]
    for i,iy in enumerate(y):
       xx.append(bin_edges[i])
       xx.append(bin_edges[i+1])
       yy.append(iy)
       yy.append(iy)

    return numpy.array(xx),numpy.array(yy)

def mf_adaptive(parts):

    masses = parts.mass.value_in(units.MSun)
    
    ones=numpy.ones_like(masses)
    m,dens=radial_density(masses, ones,N=30,dim=1)

    return m,dens

def mf_plot(parts):
#    x,y=mf_histo(parts)

    m,mf=mf_adaptive(parts)

    
    f=pyplot.figure(figsize=(8,6))
    
#    pyplot.loglog(x,y+1.e-11,'r')
    pyplot.loglog(m,mf,'b')
             
    c = ((0.1**-1.35) - (125.0**-1.35)) / 1.35
    pyplot.plot(x, len(parts)/ c * (x**-2.35))
    
    pyplot.xlim(0.1,100)
    pyplot.ylim(0.1,1.e4)
    
    fig=StringIO()
    
    pyplot.savefig(fig,format="png")
    f.clear()
    pyplot.close(f)
    return fig
  

if __name__=="__main__":

    f=pyplot.figure(figsize=(8,6))

      
    directory="./snapshots"
    i=10
      
    try:
        dm=read_set_from_file(directory+"/sink-%6.6i"%i,"amuse")
        print dm.collection_attributes.time.in_(units.Myr)
        print len(dm),dm.total_mass().in_(units.MSun),dm.mass.max().in_(units.MSun)
        nstar=len(numpy.where(dm.mass>0.1|units.MSun)[0])
        print nstar
        m,mf=mf_histo(dm)
        
        pyplot.loglog(m,mf+1.e-11,label=str(i))
        
        c = ((0.1**-1.35) - (125.0**-1.35)) / 1.35
        pyplot.plot(m, 2*nstar/ c * (m**-2.35),'r:')

    except:# Exception as ex:
        raise 
        pass
    pyplot.legend(loc='upper right', prop=dict(size=12))    
    
    pyplot.xlim(0.1,100)
    pyplot.ylim(0.1,1.e5)
    
    pyplot.show()
#    pyplot.savefig(fig,format="png")
