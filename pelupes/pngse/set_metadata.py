import pickle

from PIL import Image

from amuse.units import units

Ndata=256

img=Image.open("data%i.png"%Ndata)

from PIL import PngImagePlugin

info = PngImagePlugin.PngInfo()

f=open("minmax%i.pkl"%Ndata,"rb")
dtmin,dtmax,Mmin,Mmax,Lmin,Lmax,Tmin,Tmax=pickle.load(f, encoding="latin1")
f.close()

comment=f"""dt_min, dt_max [Myr]: {dtmin.value_in(units.Myr):.6} {dtmax.value_in(units.Myr):.6}
M_min, M_max [MSun]: {Mmin.value_in(units.MSun):.6} {Mmax.value_in(units.MSun):.6}
L_min, L_max [LSun]: {Lmin.value_in(units.LSun):.6} {Lmax.value_in(units.LSun):.6}
T_min, T_max [K]: {Tmin.value_in(units.K):.6} {Tmax.value_in(units.K):.6}
"""

info.add_text("comment",comment)

img.save("SE-%i.png"%Ndata, "png", pnginfo=info)
nimg = Image.open("SE-%i.png"%Ndata)

print(nimg.info['comment'])
