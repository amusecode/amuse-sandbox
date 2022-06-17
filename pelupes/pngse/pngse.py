# to get the png: se.py ->se_data_to_png.py -> set_metadata.py
import numpy
from PIL import Image

from amuse.units import units

from amuse.units.quantities import VectorQuantity

from matplotlib import pyplot

from amuse.datamodel import Particle, Particles

from amuse.community.sse.interface import SSE

def convert_from_byte(x,mn,mx):
    return mn*10**(x/255.*numpy.log10(mx/mn))

NMassBase=256 
Mmin=0.1 | units.MSun
Mmax=100. | units.MSun

class PNGSE(object):

    def __init__(self,Ndata=256):
        self.Ndata=Ndata
        self.ratio=NMassBase/Ndata

        self._data=Image.open("SE-%i.png"%Ndata)

        # extract limits from metadata
        comment=self._data.info['comment']
        
        comment=comment.split("\n")
        l=comment[0].split(" ")
        dtmin, dtmax=float(l[-2]) | units.Myr,float(l[-1]) | units.Myr
        l=comment[1].split(" ")
        Mmin, Mmax=float(l[-2]) | units.MSun,float(l[-1]) | units.MSun
        l=comment[2].split(" ")
        Lmin, Lmax=float(l[-2])| units.LSun ,float(l[-1]) | units.LSun
        l=comment[3].split(" ")
        Tmin, Tmax=float(l[-2]) | units.K,float(l[-1]) | units.K

        #~ print(dtmin,dtmax)
        #~ print(Mmin,Mmax)
        #~ print(Lmin,Lmax)
        #~ print(Tmin,Tmax)
        self.dtmin=dtmin
        self.dtmax=dtmax
        self.Mmin=Mmin
        self.Mmax=Mmax
        self.Lmin=Lmin
        self.Lmax=Lmax
        
        size=self._data.size

        assert size[0]==Ndata

        rdata=[]

        for mi in range(size[0]):
            age=numpy.zeros(size[1]) | units.Myr
            timestep=numpy.zeros(size[1]) | units.Myr
            mass=numpy.zeros(size[1]) | units.MSun
            lum=numpy.zeros(size[1]) | units.LSun
            temp=numpy.zeros(size[1]) | units.K
          
            tnew=0| units.Myr
            for ti in range(size[1]):
                age[ti]=tnew
                l,t,dt,m=self._data.getpixel((mi,ti))
                dt=convert_from_byte(dt,dtmin,dtmax)
                m=convert_from_byte(m,Mmin,Mmax)
                l=convert_from_byte(l,Lmin,Lmax)
                t=convert_from_byte(t,Tmin,Tmax)
                mass[ti]=m
                lum[ti]=l
                temp[ti]=t
                timestep[ti]=dt
                tnew=tnew+dt
          
            rdata.append([age,mass,lum,temp,timestep])
    
        #~ for i,(a,m,l,t) in enumerate(rdata):
          #~ print(">",i, a[0],m[0],l[0],t[0])
    
        self._rdata=rdata
    
        self.particles=Particles()
        self.model_time=0.| units.Myr

    def commit_particles(self):
        self.particles.zams_mass=self.particles.mass
        self.evolve_model(self.model_time)
    
# evolve:
# - bepaal mass bin
# - mass schaal? 
# - bepaal L schaal
# - tschaal

    def evolve_particle(self,mass,tend):
        
        mass=numpy.clip(mass.value_in(units.MSun),self.Mmin.value_in(units.MSun),
                        self.Mmax.value_in(units.MSun)) | units.MSun
        lmass=(NMassBase-1)*numpy.log10(mass/self.Mmin)/numpy.log10(self.Mmax/self.Mmin)/self.ratio
        #~ imass=int(numpy.ceil(lmass))
        imass=int(numpy.round(lmass))
        dm=lmass-imass
        mfac=10**(dm*self.ratio/(NMassBase-1.)*numpy.log10(self.Mmax/self.Mmin))  # this is just ratio between mass and mass bin imass
        
        imass1=imass+numpy.sign(dm)  # imass can also be >Ndata-1 ...
        #~ print(imass,imass1,mfac,(1/mfac*mass).in_(units.MSun),"<<")
        imass1=numpy.clip(imass1,0,self.Ndata-1)

        #~ print self.Mmin*10**(imass*self.ratio/(NMassBase-1.)*numpy.log10(self.Mmax/self.Mmin))
        #~ print self.Mmin*10**(imass1*self.ratio/(NMassBase-1.)*numpy.log10(self.Mmax/self.Mmin))

        #~ print(imass,imass1)

        if imass1!=imass:
          l1=self._data.getpixel((imass,0))[0]
          l2=self._data.getpixel((imass1,0))[0]          
          dl=(l2-l1)*dm/(imass1-imass) #
          lfac=10**(dl/(255.)*numpy.log10(self.Lmax/self.Lmin))
          tfac=lfac/mfac
        else:
          lfac=1.
          tfac=1.
        
        #~ print("evolution factors:", mfac,lfac,tfac)
        age,mass,lum,temp,timestep=self._rdata[imass]
        ts=tend*tfac
        i=1
        while i<len(age)-1 and age[i]<ts:
          i+=1
        #~ m=mfac*(mass[i-1]+(mass[i]-mass[i-1])*(ts-age[i-1])/(age[i]-age[i-1]))
        #~ l=lfac*(lum[i-1]+(lum[i]-lum[i-1])*(ts-age[i-1])/(age[i]-age[i-1]))
        #~ t=temp[i-1]+(temp[i]-temp[i-1])*(ts-age[i-1])/(age[i]-age[i-1])
        f=(ts-age[i-1])/(age[i]-age[i-1])
        if f>1: 
          f=1
        m=mfac*mass[i-1]**(1-f)*mass[i]**(f)
        l=lfac*lum[i-1]**(1-f)*lum[i]**(f)
        t=temp[i-1]**(1-f)*temp[i]**(f)
        dt=tfac*timestep[i-1]**(1-f)*timestep[i]**(f)
        #~ print(">>:", f,l, lum[i-1], lum[i])
        return m,l,t,dt

    
    def evolve_model(self,tend):
        for p in self.particles:
            mass,lum,temp,dt=self.evolve_particle(p.zams_mass,tend)
            p.age=tend
            p.mass=mass
            p.luminosity=lum
            p.temperature=temp
            p.timestep=dt
        self.model_time=tend

    def stop(self):
        pass

    def reset(self):
        self.model_time=0. | units.Myr
        self.particles=Particles()

# in the absence of stellar_remnant...
def dimm_star(star):
    return (star.luminosity<= 1.e-3 | units.LSun)

           
if __name__=="__main__":
    """ simple test against SSE """

    
    def stellar_remnant_state(star):
        return (10 <= star.stellar_type.value_in(units.stellar_type))* \
            (star.stellar_type.value_in(units.stellar_type) < 16)
  

    se=PNGSE(256)
    
    m=12.3 | units.MSun
    p=Particle(zams_mass=m)
    p2=se.particles.add_particle(p)

    sse=SSE(channel_type="sockets")
    p=Particle(mass=m)
    p1=sse.particles.add_particle(p)

    print("initial models:")
    print(p1.mass.value_in(units.MSun),p1.luminosity.value_in(units.LSun), p1.temperature)
    print(p2.mass.value_in(units.MSun),p2.luminosity.value_in(units.LSun), p2.temperature)

    done = False
    p1a=[]
    p2a=[]
    p1l=[]
    p2l=[]
    p1t=[]
    p2t=[]
    p1m=[]
    p2m=[]
    while not done:
      p1a.append(p1.age.value_in(units.Myr))
      p1t.append(p1.temperature.value_in(units.K))
      p1m.append(p1.mass.value_in(units.MSun))
      p1l.append(p1.luminosity.value_in(units.LSun))
      p2a.append(p2.age.value_in(units.Myr))
      p2t.append(p2.temperature.value_in(units.K))
      p2m.append(p2.mass.value_in(units.MSun))
      p2l.append(p2.luminosity.value_in(units.LSun))

      #~ done = stellar_remnant_state(p1)
      done = dimm_star(p2)
      
      #~ timestep=p1.time_step
      timestep=.01 | units.Myr
      p1.evolve_for(timestep)
      se.evolve_model(p1.age)
      #~ print(p1.age)
    
    pyplot.loglog(p1a,p1l,'g')
    pyplot.loglog(p2a,p2l,'b')
    pyplot.show()
