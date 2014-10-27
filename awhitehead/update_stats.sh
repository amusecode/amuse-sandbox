#!/bin/bash

python /home/draco/alf/amuse/sandbox/awhitehead/get_lagrangian_radii.py /home/draco/alf/data/output/32k_plummer_05pct.log >/home/draco/alf/data/output/32k_plummer_05pct_lagrad.csv
python /home/draco/alf/amuse/sandbox/awhitehead/plot_lagrad.py /home/draco/alf/data/output/32k_plummer_05pct_lagrad.csv

python /home/draco/alf/amuse/sandbox/awhitehead/get_lagrangian_radii.py /home/draco/alf/data/output/32k_plummer_10pct.log >/home/draco/alf/data/output/32k_plummer_10pct_lagrad.csv
python /home/draco/alf/amuse/sandbox/awhitehead/plot_lagrad.py /home/draco/alf/data/output/32k_plummer_10pct_lagrad.csv

python /home/draco/alf/amuse/sandbox/awhitehead/get_lagrangian_radii.py /home/draco/alf/data/output/32k_plummer_20pct.log >/home/draco/alf/data/output/32k_plummer_20pct_lagrad.csv
python /home/draco/alf/amuse/sandbox/awhitehead/plot_lagrad.py /home/draco/alf/data/output/32k_plummer_20pct_lagrad.csv

python /home/draco/alf/amuse/sandbox/awhitehead/get_lagrangian_radii.py /home/draco/alf/data/output/32k_plummer_30pct.log >/home/draco/alf/data/output/32k_plummer_30pct_lagrad.csv
python /home/draco/alf/amuse/sandbox/awhitehead/plot_lagrad.py /home/draco/alf/data/output/32k_plummer_30pct_lagrad.csv

python /home/draco/alf/amuse/sandbox/awhitehead/get_lagrangian_radii.py /home/draco/alf/data/output/4k_plummer_50pct.log >/home/draco/alf/data/output/4k_plummer_50pct_lagrad.csv
python /home/draco/alf/amuse/sandbox/awhitehead/plot_lagrad.py /home/draco/alf/data/output/4k_plummer_50pct_lagrad.csv

cd /home/draco/alf/data/output
cadaver -t <<EOF
open https://dav.box.com/dav
cd astro
mput 32k_plummer_*_lagrad.pdf
mput 4k_plummer_*_lagrad.pdf
quit
EOF

for QTY in NPU Ntot mass Etot Etop Eext Eint Eerr Edel Ecor Rvir Qvir Rcore Mcore Nmul Nbin Emul
do
	for PCT in '05' '10' '20' '30'
	do
		python /home/draco/alf/amuse/sandbox/awhitehead/plotq.py -f /home/draco/alf/data/output/32k_plummer_${PCT}pct.log -q $QTY
	done
	for PCT in '50'
	do
		python /home/draco/alf/amuse/sandbox/awhitehead/plotq.py -f /home/draco/alf/data/output/4k_plummer_${PCT}pct.log -q $QTY
	done
	cd /home/draco/alf/data/output
	cadaver -t <<EOF
open https://dav.box.com/dav
cd astro
mput 32k_plummer_*_${QTY}.pdf
mput 4k_plummer_*_${QTY}.pdf
quit
EOF
done

