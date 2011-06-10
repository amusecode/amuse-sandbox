clustergas scripts

boxedfi.py: not used, wrapper class for Fi to remove far away particles 
clustergas_gadget.py: actual starting routine for clustergas simulations, sets IC, makes grav_gas_sse system and evolves
clustergas.py: old script (with fi)
copycat.py: copycat class for gravity bridging
fast.py: "bridge" class with SPH
grav_gas_sse.py: bridge+sse+feedback class
lmech.py: mechanical luminosity routines (these may no longer be necessary since SE has Mdot and Vinf routines?
prof.py: parsing profiling output example
readme.txt: this readme
run.py: driver script for experiments (note: not all parameters are settable from here, some - less frequently changed-
        are hardcoded in clustergas_gadget (or even in grav_gas_sse)
SSEplus.py: SE with tracking of E_feedback



