#!/bin/bash

python /home/draco/alf/amuse/sandbox/awhitehead/get_lagrangian_radii.py /home/draco/alf/data/output/32k_plummer_05pct.log >/home/draco/alf/data/output/32k_plummer_05pct_lagrad.csv
python /home/draco/alf/amuse/sandbox/awhitehead/plot_lagrad.py /home/draco/alf/data/output/32k_plummer_05pct_lagrad.csv

python /home/draco/alf/amuse/sandbox/awhitehead/get_lagrangian_radii.py /home/draco/alf/data/output/32k_plummer_10pct.log >/home/draco/alf/data/output/32k_plummer_10pct_lagrad.csv
python /home/draco/alf/amuse/sandbox/awhitehead/plot_lagrad.py /home/draco/alf/data/output/32k_plummer_10pct_lagrad.csv

python /home/draco/alf/amuse/sandbox/awhitehead/get_lagrangian_radii.py /home/draco/alf/data/output/32k_plummer_20pct.log >/home/draco/alf/data/output/32k_plummer_20pct_lagrad.csv
python /home/draco/alf/amuse/sandbox/awhitehead/plot_lagrad.py /home/draco/alf/data/output/32k_plummer_20pct_lagrad.csv

python /home/draco/alf/amuse/sandbox/awhitehead/get_lagrangian_radii.py /home/draco/alf/data/output/32k_plummer_30pct.log >/home/draco/alf/data/output/32k_plummer_30pct_lagrad.csv
python /home/draco/alf/amuse/sandbox/awhitehead/plot_lagrad.py /home/draco/alf/data/output/32k_plummer_30pct_lagrad.csv

cd /home/draco/alf/data/output
cadaver -t <<EOF
open https://dav.box.com/dav
cd astro
mput 32k_plummer_*_lagrad.pdf
quit
EOF

