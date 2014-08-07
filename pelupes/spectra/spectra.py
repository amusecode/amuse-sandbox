from amuse.units import units,constants
import numpy
import uvblue
import bluered
from os import path

#uvdir="/disks/koppoel1/farias/spectral_data/data/uvblue/"
# take care with this, I dont have write permisions in the uvdir
# so available.txt is in brdir 
uvdir="/home/inti/code/uvbluered/UVBLUE/"
brdir="/home/inti/code/uvbluered/BLUERED/"
#brdir="/home/inti/code/uvbluered/bluered_r10/"


pi=numpy.pi
e=numpy.e
kB=constants.kB
h=constants.h
c=constants.c
Ry=constants.Rydberg_constant
sigma=constants.Stefan_hyphen_Boltzmann_constant

def interp1d(x,y):
  def f(x_):
    return numpy.interp(x_,x.copy(),y.copy())
  return f

def B_nu(nu,t):
  return 2.*pi*h*nu**3/c**2 * 1./ (e**(h*nu/kB/t)-1)

def B_lambda(l,t):
  return (2.*pi*h*c**2)/(l**5)*(1./(e**(h*c/(l*kB*t))-1))

def closest_spectrum(Teff,metallicity=0.0,eta=2,logg=5.0,uvmatch=True,brmatch=True,verbose=False):

    brOK=numpy.loadtxt(brdir+"br_available.txt")
    uvOK=numpy.loadtxt(uvdir+"uv_available.txt")

    d_uv  = ((uvOK[:,0]-Teff.value_in(units.K))**2 + (uvOK[:,1]-logg*10)**2 +(uvOK[:,2]-metallicity*10)**2)**0.5
    d_br  = ((brOK[:,0]-Teff.value_in(units.K))**2 + (brOK[:,1]-logg*10)**2 +(brOK[:,2]-metallicity*10)**2)**0.5

    if (uvmatch and brmatch):
            temp_br=abs(d_br - min(d_uv))
            temp_uv=abs(d_uv - min(d_br))
            index_br=temp_br.argmin()
            index_uv=temp_uv.argmin()

            while (brOK[index_br] != uvOK[index_uv]).any() :
                    if d_uv[index_uv] > d_br[index_br]:
                            temp_br[index_br]=max(temp_br)
                            index_br=temp_br.argmin() 
                    elif d_uv[index_uv] < d_br[index_br]:
                            temp_uv[index_uv]=max(temp_uv)
                            index_uv=temp_uv.argmin()


            uv_flag=True
            br_flag=True

    elif uvmatch:
            index_uv=d_uv.argmin()
            uv_flag=True
            br_flag=False
    elif brmatch:
            index_br=d_br.argmin()
            uv_flag=False
            br_flag=True
    else:
            if min(d_uv) < min(d_br):
                    index_uv=d_uv.argmin()
                    uv_flag=True
                    br_flag=False
            elif min(d_br) < min (d_uv):       
                    index_br=d_br.argmin()
                    uv_flag-False
                    br_flag=True
            elif min(d_br)==min(d_uv):
                    temp_br=abs(d_br - min(d_uv))
                    temp_uv=abs(d_uv - min(d_br))
                    index_br=temp_br.argmin()
                    index_uv=temp_uv.argmin()
                    while (brOK[index_br] != uvOK[index_uv]).any() :
                            if d_uv[index_uv] > d_br[index_br]:
                                    temp_br[index_br]=max(temp_br)
                                    index_br=temp_br.argmin() 
                            elif d_uv[index_uv] < d_br[index_br]:
                                    temp_uv[index_uv]=max(temp_uv)
                                    index_uv=temp_uv.argmin()
                    uv_flag=True
                    br_flag=True

    if (uv_flag and br_flag ):
            T_cl=brOK[index_br,0]
            logg_cl=brOK[index_br,1]*0.1
            met_cl=brOK[index_br,2]*0.1
            eta_cl=brOK[index_br,3]
            if verbose:
                    print "booth libraries used"
    elif (uv_flag and not br_flag):
            T_cl=uvOK[index_uv,0]
            logg_cl=uvOK[index_uv,1]*0.1
            met_cl=uvOK[index_uv,2]*0.1
            eta_cl=uvOK[index_uv,3]
            if verbose:
                    print "just UVBLUE library used"
    elif (br_flag and not uv_flag):
            T_cl=brOK[index_br,0]
            logg_cl=brOK[index_br,1]*0.1
            met_cl=brOK[index_br,2]*0.1
            eta_cl=brOK[index_br,3]
            if verbose:
                    print "just BLUERED library used"
    else:
            raise Exception("BUG: Not supported case")

    return T_cl,logg_cl,met_cl,eta_cl
    
    
def read_spectrum(Teff=50000,metallicity=0.0,eta=2,logg=5,verbose=False):

    #looking for the right file

    if metallicity >= 0:
            m="p"+'{:0>2}'.format(str(int(10*metallicity)))
    else:
            m="m"+'{:0>2}'.format(str(int(abs(10*metallicity))))
    
    T="t"+'{:0>5}'.format(str(int(Teff)))

    g="g"+'{:0>2}'.format(str(int(10*logg)))

    k="k"+str(int(eta))
    #uvfile= uvdir + m + "/" + T + g + m + k + ".flx.bz2"
    uvfile= uvdir + T + g + m + k + ".flx.gz"
    brfile= brdir + "spectrum." + T + g + m + k + ".gz"

    if not path.isfile(uvfile):
            uvfile=None
            if verbose:
                    print "UVBLUE file does not exist"
    if not path.isfile(brfile):
            brfile=None
            if verbose:
                    print "BLUERED file does not exist"
    if ((uvfile is None) and (brfile is None) ):
            raise Exception("No library file to read")

    return uvfile,brfile


def obtain_available_data():

    brOK=numpy.loadtxt(brdir+"br_available.txt")
    uvOK=numpy.loadtxt(brdir+"uv_available.txt")

    return brOK,uvOK

def obtain_available_metallicity():

     BR,UV=obtain_available_data()
     available_br=list()
     available_uv=list()
     _br=list()
     _uv=list()
     for i in range(0,len(BR)):
             if not (BR[i,2] in available_br):
                     available_br.append(BR[i,2])
                     _br.append(round(BR[i,2]*0.1,2))

     for i in range(0,len(UV)):
             if not (UV[i,2] in available_uv):
                     available_uv.append(UV[i,2])
                     _uv.append(round(UV[i,2]*0.1,2))

     _br.sort()
     _uv.sort()
     return _br,_uv

def check_metallicity(meta):

     _br,_uv=obtain_available_metallicity()
     out="ok"
     if not (meta in _br) | (meta in _uv):
             raise Exception("Not available metallicity")
     if not (meta in _br):
             print "WARNING: metallicity available just in the BLUERED database"
             out="uv"
     if not (meta in _uv):
             print "WARNING: metallicity available just in the UVBLUE database"
             out="br"
 
     return out


def plot_available_data(metallicity=0.0):

     from matplotlib import pyplot
 
     _br,_uv=obtain_available_metallicity()
     print "available metallicities:"
     print "BLUERED",_br
     print "UVBLUE", _uv
 
     meta=metallicity*10.
 
     check=check_metallicity(metallicity)
     BR,UV=obtain_available_data()
 
     br=list()
     uv=list()
     for i in range(0,len(BR)):
             if BR[i,2]==meta :
                     br.append(BR[i,:])
 
     for i in range(0,len(UV)):
             if UV[i,2]==meta :
                     uv.append(UV[i,:])
 
     
     br=numpy.array(br)
     uv=numpy.array(uv)

     if (check=="br"):
             uv=br*0.
     elif (check=="uv"):
             br=uv*0.

     pyplot.xlabel("Temperature [K]")
     pyplot.ylabel("log(g)")

     Tmax=max(max(uv[:,0]),max(br[:,0]))
     Tmin=min(min(uv[:,0]),min(br[:,0])) 

     gmax=0.1*max((max(uv[:,1]),max(br[:,1])))
     gmin=0.1*min((min(uv[:,1]),min(br[:,1])))
 
     margin=0.1
 
     pyplot.xlim((Tmax+margin*Tmax,Tmin-margin*Tmax))
     pyplot.ylim((gmax+margin*gmax,gmin-margin*gmax))
 
     pyplot.text(Tmax*0.5,gmax*0.0,"metallicity=%2.1f" % round(meta*0.1,2),ha="center")
     if (check=="ok"):
             pyplot.plot(uv[:,0],uv[:,1]/10.,"ob",label="UVBLUE")
             pyplot.plot(br[:,0],br[:,1]/10.,"or",label="BLUERED")
     elif (check=="br"):
             pyplot.plot(br[:,0],br[:,1]/10.,"or",label="BLUERED")
     elif (check=="uv"):
             pyplot.plot(uv[:,0],uv[:,1]/10.,"ob",label="UVBLUE")

     pyplot.legend(loc=2)

class Spectrum:
        def __init__(self,star,
                        metallicity=0.0,
                        eta=2,
                        uvmatch=True,
                        brmatch=True,
                        lambda_resolution=None,
                        resolution=None,
                        verbose=False):

            self.Teff=star.temperature
            self.logg=numpy.log10(star.mass.value_in(units.MSun)) - 2.*numpy.log10(star.radius.value_in(units.RSun)) + 4.437
            self.metallicity=metallicity
            self.eta=eta

            self.get_spectrum(uvmatch=uvmatch,
                            brmatch=brmatch,
                            resolution=resolution,
                            lambda_resolution=lambda_resolution,
                            verbose=verbose)


        def get_spectrum(self,uvmatch=True,brmatch=True,resolution=None,lambda_resolution=None,verbose=False):

            T_ref,lg_ref,met_ref,eta_ref=closest_spectrum(Teff=self.Teff,
                                         metallicity=self.metallicity,
                                         eta=self.eta,
                                         logg=self.logg,
                                         uvmatch=uvmatch,
                                         brmatch=brmatch)
            if verbose:
                    print "Closest spectrum in the library:"
                    print "Teff  : %5.2f [K]" %T_ref
                    print "metallicity : %5.2f" %met_ref
                    print "log(g) : %5.2f" %lg_ref
                    print "eta : %5.2f" %eta_ref
                    print ""

            self.Tref=units.K(T_ref)
            self.logg_ref=lg_ref
            self.metallicity_ref=met_ref
            self.lambda_resolution=lambda_resolution
            self.resolution=resolution

            uvfile,brfile=read_spectrum(Teff=T_ref,
                                        metallicity=met_ref,
                                        eta=eta_ref,
                                        logg=lg_ref,
                                        verbose=verbose)
            
            uvstart=850.
            brend=7000.


            if ((uvfile) and (brfile)):
                    brdata=bluered.UVBlue(brfile)
                    uvdata=uvblue.UVBlue(uvfile)

                    res=min(uvdata.resolution,brdata.resolution)

                    if resolution==None:
                            self.resolution=res
                    elif resolution > res:
                            if verbose:
                                    print "Warning: input resolution had to be smaller than", res
                                    print "Using", res,"instead"
                            self.resolution=res

                    elif resolution <= res:
                            self.resolution=resolution

                    lres_uv= round( uvdata.lamb[1] / (uvdata.lamb[1]-uvdata.lamb[0]),ndigits=1)
                    lres_br= round( brdata.lamb[1] / (brdata.lamb[1]-brdata.lamb[0]),ndigits=1)

                    l_res=min(lres_uv,lres_br)

                    if lambda_resolution==None:
                            self.lambda_resolution=l_res
                    elif lambda_resolution > l_res:
                            if verbose:
                                    print "Warning: input lambda resolution had to be smaller than", l_res
                                    print "Using", l_res,"instead"
                            self.lambda_resolution=l_res

                    elif lambda_resolution <= l_res:
                            self.lambda_resolution=lambda_resolution

                    brstart=min(brdata.lamb)
                    uvend=max(uvdata.lamb)

                    end_index=min(min(numpy.where(brdata.lamb>=uvend)))

                    if lres_uv > self.lambda_resolution :
                            li=brdata.lamb[end_index]
                            l_uv=list()
                            while li > (uvstart + li/self.lambda_resolution):
                                    li-=li/self.lambda_resolution
                                    l_uv.append(li)
                            l_uv=numpy.sort(numpy.array(l_uv))

                    if lres_br > self.lambda_resolution :
                            st_index=min(min(numpy.where(l_uv>=brstart)))
                            li=l_uv[st_index]
                            l_br=list()
                            l_br.append(li)
                            while li <= (brend - li/self.lambda_resolution):
                                    li+=li/self.lambda_resolution
                                    l_br.append(li)
                            l_br=numpy.array(l_br)
                    else:
                            l_br=brdata.lamb


                    if self.resolution < uvdata.resolution:
                      if verbose:
                              print "broaden UVBLUE spectrum from R= %5.2f to R= %5.2f" %(uvdata.resolution,self.resolution)
                      uv_flam=uvblue.broaden_spectrum(uvdata.lamb,uvdata.fnlam,res=self.resolution,_res=uvdata.resolution)
                    else:
                      uv_flam=uvdata.fnlam

                    uv_inter=interp1d(uvdata.lamb,uv_flam)

                    
                    if self.resolution < brdata.resolution:
                            if verbose:
                                    print "broaden BLUERED spectrum from R= %5.2f to R= %5.2f" %(brdata.resolution,self.resolution)
                            br_flam=bluered.broaden_spectrum(brdata.lamb,brdata.fnlam,res=self.resolution,_res=brdata.resolution)
                    else:
                            br_flam=brdata.fnlam

                    br_inter=interp1d(brdata.lamb,br_flam)


                    br_flam=br_inter(l_br)
                    uv_flam=uv_inter(l_uv)

                    self.uv_flam=uv_flam
                    self.uv_inter=uv_inter
                    self.uv_lamb=l_uv

                    ######## combine spectra

                    start_index=min(min(numpy.where(l_uv >= brstart))) - 2
                    end_index=max(max(numpy.where(l_br <= uvend))) - 1
                    lamb_shared=l_br[:end_index]
                    
                    uvend=lamb_shared[-1]

                    fr_left=abs(uvend-lamb_shared)/(uvend-brstart) 
                    fr_right=abs(brstart-lamb_shared)/(uvend-brstart)

                    flam_shared=(br_inter(lamb_shared)*fr_right + uv_inter(lamb_shared)*fr_left)
                    flam=numpy.concatenate( (uv_inter(l_uv[:start_index+1]) , flam_shared , br_inter(l_br[end_index:])) )
                    flam=flam*5.6697e-5*(self.Tref.value_in(units.K)**4) 
                    lamb=numpy.concatenate((l_uv[:start_index+1],l_br)) | units.angstrom
                    #########################

                    #scaling
                    BB_teff=B_lambda(lamb,self.Teff)
                    BB_tref=B_lambda(lamb,self.Tref)
                    scale=BB_teff/BB_tref
                    flam=flam*scale | units.erg*units.s**-1*units.cm**-2*units.angstrom**-1

                    self.uv_flag=True
                    self.br_flag=True

            elif (uvfile and (brfile==None)):
                    int_flag=False
                    uvdata=uvblue.UVBlue(uvfile)
                    uvend=max(uvdata.lamb)
                    uvstart=min(uvdata.lamb)

                    if lambda_resolution==None:
                            self.lambda_resolution= round( uvdata.lamb[1] / (uvdata.lamb[1]-uvdata.lamb[0]),ndigits=1)
                            l_uv=uvblue.lamb
                    else:
                            self.lambda_resolution=lambda_resolution

                            l_uv=list()
                            li=uvend
                            l_uv.append(li)
                            while li > uvstart+li/self.lambda_resolution:
                                    li-=li/self.lambda_resolution
                                    l_uv.append(li)
                            l_uv=numpy.sort(numpy.array(l_uv))
                            int_flag=True


                    if resolution==None or resolution==uvdata.resolution:
                            uv_flam=uvdata.fnlam
                            self.resolution=uvdata.resolution
                    elif resolution <= uvdata.resolution:
                            self.resolution=resolution
                            if verbose:
                                    print "broaden UVBLUE spectrum from R= %5.2f to R= %5.2f" %(uvdata.resolution,self.resolution)
                            uv_flam=uvblue.broaden_spectrum(uvdata.lamb,uvdata.fnlam,res=self.resolution)
                            int_flag=True

                    if int_flag:
                            uv_inter=interp1d(uvdata.lamb,uv_flam)
                            flam_uv=uv_inter(l_uv)
                    else:
                            flam_uv=uvdata.fnlam
                    


                    l_br=list()
                    li=uvend

                    while li <= brend:
                            li+=li/self.lambda_resolution
                            l_br.append(li)
                    l_br=numpy.array(l_br)

                    lamb=numpy.concatenate((l_uv,l_br)) |units.angstrom

                    BB_teff=B_lambda(units.angstrom(l_uv),self.Teff)
                    BB_tref=B_lambda(units.angstrom(l_uv),self.Tref)
                    scale=BB_teff/BB_tref

                    flam_br=B_lambda(units.angstrom(l_br),self.Teff).value_in(units.erg*units.s**-1*units.cm**-2*units.angstrom**-1)
                    flam_uv=flam_uv*5.6697e-5*(self.Tref.value_in(units.K)**4) 

                    flam=numpy.concatenate((flam_uv*scale,flam_br)) | units.erg*units.s**-1*units.cm**-2*units.angstrom**-1

                    self.uv_flag=True
                    self.br_flag=False

            elif (brfile and (uvfile==None)):
                    int_flag=False
                    brdata=bluered.UVBlue(brfile)
                    brend=max(brdata.lamb)
                    brstart=min(brdata.lamb)

                    if lambda_resolution==None:
                            self.lambda_resolution= round( brdata.lamb[1] / (brdata.lamb[1]-brdata.lamb[0]),ndigits=1)
                            l_br=brblue.lamb
                    else:
                            self.lambda_resolution=lambda_resolution

                            l_br=list()
                            li=brstart
                            l_br.append(li)
                            while li < brend-li/self.lambda_resolution:
                                    li+=li/self.lambda_resolution
                                    l_br.append(li)
                            l_br=numpy.array(l_br)
                            int_flag=True


                    if resolution==None or resolution==brdata.resolution:
                            br_flam=brdata.fnlam
                            self.resolution=brdata.resolution
                    elif resolution <= brdata.resolution:
                            self.resolution=resolution
                            if verbose:
                                    print "broaden BLUERED spectrum from R= %5.2f to R= %5.2f" %(brdata.resolution,self.resolution)
                            br_flam=bluered.broaden_spectrum(brdata.lamb,brdata.fnlam,res=self.resolution)
                            int_flag=True

                    if int_flag:
                            br_inter=interp1d(brdata.lamb,br_flam)
                            flam_br=br_inter(l_br)
                    else:
                            flam_br=brdata.fnlam
                    


                    l_uv=list()
                    li=brstart

                    while li >= uvstart:
                            li-=li/self.lambda_resolution
                            l_uv.append(li)
                    l_uv=numpy.sort(numpy.array(l_uv))

                    lamb=numpy.concatenate((l_uv,l_br)) |units.angstrom

                    BB_teff=B_lambda(units.angstrom(l_br),self.Teff)
                    BB_tref=B_lambda(units.angstrom(l_br),self.Tref)
                    scale=BB_teff/BB_tref

                    flam_uv=B_lambda(units.angstrom(l_uv),self.Teff).value_in(units.erg*units.s**-1*units.cm**-2*units.angstrom**-1)
                    flam_br=flam_br*5.6697e-5*(self.Tref.value_in(units.K)**4) 

                    flam=numpy.concatenate((flam_uv,flam_br*scale)) | units.erg*units.s**-1*units.cm**-2*units.angstrom**-1

                    self.uv_flag=False
                    self.br_flag=True
            else:
                    raise Exception("no data selected")
             
            self.flux=flam 
            self.lamb=lamb.in_(units.angstrom)
