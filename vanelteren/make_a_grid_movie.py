import sys
import unittest
import numpy 
import random
import collections
import os
import subprocess
import time
import gc


try:
    from matplotlib import pyplot
    from matplotlib import colors
    from matplotlib import colorbar
    from matplotlib import cm
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from amuse import io
from amuse.units import generic_unit_system
from amuse.units import units

class PlotGridRho(object):
    
    def __init__(self):
        self.length_units = generic_unit_system.length
        self.mass_units = generic_unit_system.mass
        self.density_units = generic_unit_system.density
        self.cmap = None
        self.norm = None
        
    def plot(self, plot, grid):
        x_values = grid.x.value_in(self.length_units)
        y_values = grid.y.value_in(self.length_units)
        rho_values = grid.rho.value_in(self.density_units)
        if self.norm is None:
            
            self.norm = colors.LogNorm(vmin=rho_values.min(), vmax=rho_values.max())
            self.cmap = cm.get_cmap('jet', 1024)
            
        im = plot.imshow(
            rho_values, 
            origin = 'lower',
            cmap=self.cmap, 
            norm=self.norm
        )
        colorbar.Colorbar(plot, im)
        

        
        
        
class MakeAMovie(object):
    
    
    def __init__(self, name_of_the_hdf_file):
        self.grid = io.read_set_from_file(name_of_the_hdf_file, 'amuse')
        
    def start(self):
        
        index = 0
        
        make_plot_objects = []
        make_plot_objects.append(PlotGridRho())
        
        number_of_plots = len(list(self.grid.history))
        t0 = time.time()
        
        print "starting..."
        figure = pyplot.figure(figsize=(8,2))
        for data in self.grid.history:
        
            for i, make_a_plot in enumerate(make_plot_objects):
                plot = figure.add_subplot(1,1,i+1)
                make_a_plot.plot(plot, data)
                
            index += 1
            figure.savefig("frame_"+format(index, "05d")+".png")
            t1 = time.time()
            dt = t1 - t0
            estimated_number_of_seconds_to_go = (dt / (index + 1)) * (number_of_plots-index)
            print index, "/", number_of_plots, " estimate to be finished in ", estimated_number_of_seconds_to_go, " seconds"
            
            figure.clear()
        
        pyplot.close(figure)
            
        #for i in range(number_of_plots):
        #    filename1 = "frame_"+format(i+1, "05d")+".png"
        #    filename2 = "frame_"+format(i+1, "05d")+".jpg"
        #    subprocess.call([
        #        "convert",
        #        filename1,
        #        filename2,
        #    ])
        #    os.remove(filename1)
            
            
        with open("frames.txt", "w") as f:
            for i in range(number_of_plots):
                filename = "frame_"+format(i+1, "05d")+".png"
                f.write(filename + "\n")
        
        
        subprocess.call( [
            "mencoder",
            "mf://@frames.txt",
            "-mf",
            "w=480:h=480:fps=10:type=png",
            "-ovc",
            "lavc",
            "-lavcopts",
            "vcodec=mpeg4:mbd=2:trell",
            "-oac",
            "copy", 
            "-o",
            "movie-cloudshock-10fps.avi"
        ])
        
        
        subprocess.call( [
            "mencoder",
            "mf://@frames.txt",
            "-mf",
            "w=480:h=480:fps=25:type=png",
            "-ovc",
            "lavc",
            "-lavcopts",
            "vcodec=mpeg4:mbd=2:trell",
            "-oac",
            "copy", 
            "-o",
            "movie-cloudshock-25fps.avi"
        ])
            
        
        for i in range(number_of_plots):
            filename = "frame_"+format(i+1, "05d")+".png"
            os.remove(filename)
            
if __name__ == '__main__':
    use_case = MakeAMovie(sys.argv[1])
    use_case.start()
