function h,z,omegam,omegal

 return, sqrt(omegal+(1-omegal-omegam)*(1.+z)^2+omegam*(1.+z)^3)
 
end

pro tophat, N=N,M=M, z0=z0, zcollaps=zcollaps,omegam=omegam,$
            omegab=omegab,H0=H0,RA=RA,gasN=gasN,lambda=lambda,$
	    iscomov=iscomov,center_on_max=center_on_max,dd=dd,ssfac=ssfac, $
	    readdata=readdata,writedata=writedata

; code to make tophat halo model
;
; N= basegrid for DM particles, this determines the number of particles.    
; M= mass of the halo, 
; z0= redshift for the initial conditions,
; zcollaps= target collapse redshift, 
; omegam,omegab,H0= cosmology parameters (omegam+omagelambda=1 assumed)
; RA= for output of radius (optional), 
; gasN= basegrid for gas (Ngas=0 for pure DM),
; iscomov=use comoving coord? 
; center_on_max= centre the tophat on maximum of density field?
; dd= for output of density field, 
; ssfac=normalization of perturbations (defaults to 1., optional), 
; lambda=angular momentum parameter

; the following is to easily reproduce the same initial conditions at 
; different resolution:
;  readdata=readin previous fourier data, 
;  writedata=write fourier data
; so generate high resolution with writedata on, then a lower resolution
; version with readdata  to get the corresponding initial conditions at lower
; resolution

; note that the input mass and the internal calculations are done in a
; system of units with G=1, mass unit =10^9 Msun, length unit =1 kpc. 
; but the output is converted (in a few lines at the end, marked) to gadget 
; units..  

; the output is stored in a common block, which is also used in wgadget
; routine
; note that I don't include u, specify this seperately in the initial
; conditions

; example calls:
;tophat, N=20,M=0.002, z0=100., zcollaps=30,omegam=1.0,omegab=0.05,H0=50.,gasN=0
;tophat,N=116,M=0.1,z0=99.,zcollaps=16.,omegam=0.3,omegab=0.04,H0=70.,gasN=58,lambda=0.03,/iscomov
;tophat,N=192,M=10.,z0=99.,zcollaps=10.,omegam=0.3,omegab=0.04,H0=70.,gasN=96,lambda=0.03,/iscomov,/center

; example write:
; wgadget,'comovtophat_N192_M10_z99_zc10_Ng96_l03c'

; F.I. Pelupessy 18/07/2007

;-----

common headerinfo, tnow,npart,massarr,comoving
common particleinfo, mass,pos,vel

;-----

 print,'3D tophat model'
 print,'mass (msun):',M*1.e9

if(n_elements(center_on_max) eq 0) then center_on_max=0
if(n_elements(lambda) eq 0) then lambda=0.
if(n_elements(gasN) eq 0) then gasN=0
if(n_elements(iscomov) eq 0) then iscomov=0
if(n_elements(ssfac) eq 0) then ssfac=1.
if(n_elements(readdata) eq 0) then readdata=0
if(n_elements(writedata) eq 0) then writedata=0


 tscale=978/14.9     ; these define some miscallenous scalings
 mscale=10.
 uscale=104.
 utempfac=4.86e5
 vscale=65.6
 mpcscale=0.001
 power=-2.5


; let's first solve the 1d program ( this give you the exact
; overdensity (sigma) and hubble law (heff) for the tophat model
; at the initial redshift  
tophatcollaps,zcollaps,z0,omegam,1.-omegam,sigma=sigma,heff=heff

; (this is better than taking linear approximation:)
;sigma=1+1.686*(1+zcollaps)/(1+z0)

 print,'sigma_lin,sigma=',1+1.686*(1+zcollaps)/(1+z0),sigma

 N=long(N)
 gasN=long(gasN)
 fdark=(omegam-omegab)/omegam     ; fraction of dm
 if(gasN eq 0) then fdark=1.        
 
 Pi=!DPi
 H=H0/vscale*mpcscale*h(z0,omegam,1.-omegam)
 Heff=H0*Heff/vscale*mpcscale
 rhocrit=3*H^2/8./Pi               ; this is the universe critical dens
 
; print,'H=',H
; print,'Heff=',Heff
; print,'rhocrit=',rhocrit
 
 rho=sigma*rhocrit                  ; this is the actual density of the 
                                    ;  spherical overdensity
; print,'rho=',rho 
 RA=(3*M/4/Pi/rho)^.333333          ; thus the radius
; print,'R=',RA 
 R0=(3*M/4/Pi/rhocrit)^.333333       ; the original radius


; normalization according to bromm 2002
anorm=0d  
for i=-N/2+1,N/2 do begin
 for j=-N/2+1,N/2 do begin
  for k=-N/2+1,N/2 do begin
  if(i ne 0 OR i ne 0 OR k ne 0) then anorm=anorm+(k^2+i^2+j^2)^(power/2.)
  endfor
 endfor
endfor
; print,'sig02=',(1+zcollaps)^2/(1+z0)^2 
 anorm=(1+zcollaps)^2/(1+z0)^2/anorm
 anorm=anorm*ssfac
; print,'A=',anorm
 
; check=dblarr(N,N)
 dk=dcomplexarr(N,N,N)      ; density and displacement arrays
 fi=dcomplexarr(N,N,N)
 fj=dcomplexarr(N,N,N)
 fk=dcomplexarr(N,N,N)
 
seed=26711L
;seed=2611L

; fill the fourier grids, possibly using disk data

; ----
 if(readdata) then begin
  openr,1,'power'
  NN=N
  readu,1,NN
  if(NN LT N) then begin
   print,' error reading power data'
   stop
  endif  
  dklarge=dcomplexarr(NN,NN,NN)
  readu,1,dklarge
  close,1
  dk=dklarge(NN/2-N/2:NN/2-1+N/2,NN/2-N/2:NN/2-1+N/2,NN/2-N/2:NN/2-1+N/2)
 endif
;---------

for i=0L,N-1 do begin
 for j=0L,N-1 do begin
  for k=0L,N-1 do begin
   ki=(i-N/2+1)
   kj=(j-N/2+1)
   kk=(k-N/2+1)
   if(not(readdata)) then begin
   phi=randomu(seed)*2*Pi
   ak=randomu(seed)
   ak=sqrt(-4/Pi*alog(ak)*anorm*(ki^2+kj^2+kk^2+1.d-30)^(power/2.))
;   ak=sqrt(anorm*(ki^2+kj^2+kk^2+1.d-30)^(-1.5))
    dk(i,j,k)=ak*exp(dcomplex(0.,phi))
   endif
   fi(i,j,k)=dk(i,j,k)*complex(0,1)*ki/(ki^2+kj^2+kk^2+1.d-30)
   fj(i,j,k)=dk(i,j,k)*complex(0,1)*kj/(ki^2+kj^2+kk^2+1.d-30)
   fk(i,j,k)=dk(i,j,k)*complex(0,1)*kk/(ki^2+kj^2+kk^2+1.d-30)
  endfor
 endfor
endfor
; ----- 
 if (writedata) then begin
  openw,1,'power'
  writeu,1,N
  writeu,1,dk
  close,1
 endif
;----
; print,'check!', dk(N/2-1,N/2-1,N/2-1)
 dk(N/2-1,N/2-1,N/2-1)=0.
 fi(N/2-1,N/2-1,N/2-1)=0.
 fj(N/2-1,N/2-1,N/2-1)=0.
 fk(N/2-1,N/2-1,N/2-1)=0.
for i=0L,N-1 do begin
 for j=0L,N-1 do begin
  for k=0L,N-1 do begin
    ii=N-i-2
    jj=N-j-2
    kk=N-k-2 
    if(ii eq -1) then ii=N-1 
    if(jj eq -1) then jj=N-1 
    if(kk eq -1) then kk=N-1 
   if( dk(i,j,k) ne conj(dk(ii,jj,kk))) then begin
       dk(i,j,k)=conj(dk(ii,jj,kk))
       fi(i,j,k)=conj(fi(ii,jj,kk))
       fj(i,j,k)=conj(fj(ii,jj,kk))
       fk(i,j,k)=conj(fk(ii,jj,kk))
   endif    
  endfor
 endfor
endfor
for i=0,1 do begin
 for j=0,1 do begin
  for k=0,1 do begin
     ii=N/2-1+i*N/2
     jj=N/2-1+j*N/2
     kk=N/2-1+k*N/2
     dk(ii,jj,kk)=(dk(ii,jj,kk)+conj(dk(ii,jj,kk)))/2.
     fi(ii,jj,kk)=(fi(ii,jj,kk)+conj(fi(ii,jj,kk)))/2.
     fj(ii,jj,kk)=(fj(ii,jj,kk)+conj(fj(ii,jj,kk)))/2.
     fk(ii,jj,kk)=(fk(ii,jj,kk)+conj(fk(ii,jj,kk)))/2.
  endfor
 endfor
endfor

dk=shift(dk,-N/2+1,-N/2+1,-N/2+1)
fi=shift(fi,-N/2+1,-N/2+1,-N/2+1)
fj=shift(fj,-N/2+1,-N/2+1,-N/2+1)
fk=shift(fk,-N/2+1,-N/2+1,-N/2+1)

; ok, now do the ffts
x=RA/Pi*fft(fi,/inverse,/double)
mx=max(abs(imaginary(x)))
if(mx gt 1.d-10*RA) then print,'warning1!', mx
y=RA/Pi*fft(fj,/inverse,/double)
mx=max(abs(imaginary(y)))
if(mx gt 1.d-10*RA) then print,'warning2!', mx
z=RA/Pi*fft(fk,/inverse,/double)
mx=max(abs(imaginary(z)))
if(mx gt 1.d-10*RA) then print,'warning3!', mx

x=double(x)
y=double(y)
z=double(z)

ix=0
iy=0
iz=0
dd=RA/Pi*fft(dk,/inverse,/double)


; before we cut the spherical region we can choose to center on the 
; biggest density peak..
if(center_on_max) then begin
 mx=max(double(dd),imx)
 ix=imx mod N
 iz=imx/(N*N)
 iy=(imx-iz*N*N)/N
 if(writedata) then begin
  openw,1,'shift'
  printf,1,ix,iy,iz
  close,1
 endif
 if(readdata) then begin
  print,'*',ix,iy,iz
  openr,1,'shift'
  readf,1,ix,iy,iz
  close,1
  ix=(ix*N)/NN
  iy=(iy*N)/NN
  iz=(iz*N)/NN
  print,'*',ix,iy,iz
 endif 
 print,mx,dd(ix,iy,iz)
 x=shift(x,N/2-ix,N/2-iy,N/2-iz)
 y=shift(y,N/2-ix,N/2-iy,N/2-iz)
 z=shift(z,N/2-ix,N/2-iy,N/2-iz)
 print,dd(N/2,N/2,N/2)
endif
 dd=shift(double(dd),N/2-ix,N/2-iy,N/2-iz)
 tvscl,dd(*,*,N/2)

; set velocities
vx=x*H
vy=y*H
vz=z*H

; add displacements to the actual particle positions
for i=0,N-1 do begin
 for j=0,N-1 do begin
  for k=0,N-1 do begin
    x(i,j,k)=x(i,j,k)+i*2*RA/N+RA/N
    y(i,j,k)=y(i,j,k)+j*2*RA/N+RA/N
    z(i,j,k)=z(i,j,k)+k*2*RA/N+RA/N
    if (x(i,j,k) lt 0.) then x(i,j,k)=x(i,j,k)+2*RA
    if (y(i,j,k) lt 0.) then y(i,j,k)=y(i,j,k)+2*RA
    if (z(i,j,k) lt 0.) then z(i,j,k)=z(i,j,k)+2*RA
    if (x(i,j,k) gt 2*RA) then x(i,j,k)=x(i,j,k)-2*RA
    if (y(i,j,k) gt 2*RA) then y(i,j,k)=y(i,j,k)-2*RA
    if (z(i,j,k) gt 2*RA) then z(i,j,k)=z(i,j,k)-2*RA
  endfor
 endfor
endfor

; flatten the arrays
x=reform(x,N^3,/overwrite)
y=reform(y,N^3,/overwrite)
z=reform(z,N^3,/overwrite)

; select particles within spherical region
a=where((x-RA)^2+(y-RA)^2+(z-RA)^2 le RA^2)

nhalo=n_elements(a)

x=x(a)-RA
y=y(a)-RA
z=z(a)-RA

; add hubble expansion
vx=vx(a)+Heff*x
vy=vy(a)+Heff*y
vz=vz(a)+Heff*z

if (iscomov) then begin
vx=vx-H*x
vy=vy-H*y
vz=vz-H*z
endif

vx=vx-mean(vx)
vy=vy-mean(vy)
vz=vz-mean(vz)

massarr=dblarr(6)

mass=dblarr(nhalo)+fdark*m/nhalo
massh=fdark*m/nhalo
massarr(1)=massh


 print,'Nhalo=',nhalo

; do gas
ngas=0
if(gasN gt 0) then begin

 xgas=dblarr(gasN,gasN,gasN) 
 ygas=dblarr(gasN,gasN,gasN) 
 zgas=dblarr(gasN,gasN,gasN) 
 for i=0,gasN-1 do begin
 for j=0,gasN-1 do begin
  for k=0,gasN-1 do begin
    xgas(i,j,k)=i*2*RA/gasN+RA/gasN-RA
    ygas(i,j,k)=j*2*RA/gasN+RA/gasN-RA
    zgas(i,j,k)=k*2*RA/gasN+RA/gasN-RA
  endfor
 endfor
endfor

xgas=reform(xgas,gasN^3,/overwrite)
ygas=reform(ygas,gasN^3,/overwrite)
zgas=reform(zgas,gasN^3,/overwrite)

a=where((xgas)^2+(ygas)^2+(zgas)^2 le RA^2)

ngas=n_elements(a)
print,'Ngas=',ngas

xgas=xgas(a)
ygas=ygas(a)
zgas=zgas(a)

vxgas=Heff*xgas
vygas=Heff*ygas
vzgas=Heff*zgas
if (iscomov) then begin
vxgas=(Heff-H)*xgas
vygas=(Heff-H)*ygas
vzgas=(Heff-H)*zgas
endif

x=[xgas,x]
y=[ygas,y]
z=[zgas,z]
vx=[vxgas,vx]
vy=[vygas,vy]
vz=[vzgas,vz]

massg=dblarr(ngas)+(1-fdark)*m/ngas

mass=[massg,mass]
massarr(0)=0.

endif

npart=lonarr(6)
npart(0)=ngas
npart(1)=nhalo

tnow=0.


; add solid body rotation
Ek=total(.5*mass*(vx^2+vy^2+vz^2))
Ep=-3./5.*M^2/RA
Iang=2./5.*M*RA^2
;print,'Ek=',Ek
;print,'Ep=',Ep
;print,'Etot=',Ek+Ep
;print,'I=',Iang

lamb=lambda
om=lamb*M^2.5/Iang/sqrt(abs(Ep+Ek))

if( lambda gt 0) then print,'omega=',om
vx=vx+y*om
vy=vy-x*om

if (iscomov) then begin
 a=1./(1.+z0)
 vx=(vx)/sqrt(a)
 vy=(vy)/sqrt(a)
 vz=(vz)/sqrt(a)
 x=x/a*H0/100.
 y=y/a*H0/100.
 z=z/a*H0/100.
 redshift=double(z0)
 tnow=1./(redshift+1)
 mass=mass*H0/100.
 massarr=massarr*H0/100.
endif

 nbodies=nhalo+ngas
 pos=fltarr(3,nbodies)
 vel=fltarr(3,nbodies)
 

; this prepares the output for gadget files 
 mass=0.
 if(ngas gt 0) then mass=massg
 mass=mass/mscale
 massarr=massarr/mscale
 
 pos(0,*)=x
 pos(1,*)=y
 pos(2,*)=z
 vel(0,*)=vx*vscale
 vel(1,*)=vy*vscale
 vel(2,*)=vz*vscale
 comoving=iscomov

end

; -------------------------------------------
; routine to write out initial conditions in gadget format
pro wgadget,filename

common headerinfo, tnow,npart,massarr,comoving
common particleinfo, mass,pos,vel

IF N_params() LT 1 THEN $
  message,'Syntax - rgadget, Filename'

if(n_elements(massarr) ne 6) then message,' massarr not good'
if(n_elements(npart) ne 6) then message,' npart not good'

print, 'comoving:',comoving
  	 
time=double(tnow)
redshift=double(0.)
if(comoving) then redshift=1./tnow-1.
print,redshift 
flag_sfr=0L
flag_feedback=0L
npartTotal=npart	
bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
la=intarr(bytesleft/2)

openw,daunit,filename,/f77_unformatted,/get_lun

writeu,daunit,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la

 nbodies=npart(0)+npart(1)+npart(2)+npart(3)+npart(4)+npart(5)

 if(nbodies*3 ne n_elements(pos) or nbodies*3 ne n_Elements(vel)) then $
    message,'wrong number of pos or vel elements'

 tmp=float(pos)
 writeu,daunit,tmp 
 tmp=float(vel)
 writeu,daunit,tmp 

  origindex=lonarr(nbodies)
  for i=0L,nbodies-1 do begin
   origindex(i)=i
  endfor 

 tmp=long(origindex)
 writeu,daunit,tmp

 nwithmass=0
 if(massarr(0) eq 0 and npart(0) gt 0) then nwithmass=nwithmass+npart(0)
 if(massarr(1) eq 0 and npart(1) gt 0) then nwithmass=nwithmass+npart(1)
 if(massarr(2) eq 0 and npart(2) gt 0) then nwithmass=nwithmass+npart(2)
 if(massarr(3) eq 0 and npart(3) gt 0) then nwithmass=nwithmass+npart(3)
 if(massarr(4) eq 0 and npart(4) gt 0) then nwithmass=nwithmass+npart(4)
 if(massarr(5) eq 0 and npart(5) gt 0) then nwithmass=nwithmass+npart(5)
 
 print, nwithmass,n_elements(mass),mass
 if(nwithmass gt 0 and nwithmass ne n_elements(mass)) then message,' mass array mismatch'

 tmp=float(mass)
 if(nwithmass  gt 0) then writeu,daunit,tmp

 close,daunit
 return
end

; solver for the 1d tophat
function thsinsolve,x
 thnew=!Pi
 nmax=100
 i=0
 repeat begin
 i=i+1
 th=thnew
 thnew=th-(th-sin(th)-x)/(1-cos(th))
 endrep until abs(thnew-th) lt 1.e-5  or i gt nmax
 if(thnew-sin(thnew)-x gt 0.001) then print,'th root warning'
 return, thnew
end

; brief solution of the 1d tophat collapse
; rmax= r(tmax)/a(z_c) !!
pro tophatcollaps,zc,z0,om,oml,sigma=sigma,heff=heff
 hzc=h(zc,om,oml)
 tzc=2./3./hzc
 tz0=tzc*((1+zc)/(1+z0))^1.5
 tmax=tzc/2.
 thz0=thsinsolve(tz0/tmax*!Pi)
 rmax=4*(tmax/tzc/6./!Pi)^(2./3.)
 rz0=rmax/2.*(1-cos(thz0))
 heff=rmax/2.*sin(thz0)*!Pi/tmax/(1-cos(thz0))/rz0
 sigma=(tz0/tzc)^2/rz0^3
end 

;-------------------------------------------------------------
; main program to initialize commons

common headerinfo, tnow,npart,massarr,comoving
common particleinfo, mass,pos,vel

end

