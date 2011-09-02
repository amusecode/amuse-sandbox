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
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.io import store
from amuse.units import nbody_system
from amuse.units import units

class PlotXYandMass(object):
    
    def __init__(self, indices_to_plot = None):
        self.length_units = nbody_system.length
        self.mass_units = nbody_system.mass
        self.indices_to_plot = indices_to_plot
        
    def plot(self, plot, particles):
        x_values = particles.x.value_in(self.length_units)
        y_values = particles.y.value_in(self.length_units)
        mass_values = particles.mass.value_in(self.mass_units)
        
    
        if not self.indices_to_plot is None:
            if not self.indices_to_plot:
                plot.set_xlim(-4.0, 4.0)
                plot.set_ylim(-4.0, 4.0)
                return
                
            x_values = x_values.take(self.indices_to_plot)
            y_values = y_values.take(self.indices_to_plot)
            mass_values = mass_values.take(self.indices_to_plot)
        
        center_of_mass = particles.center_of_mass()
        
        x_values -= center_of_mass.x.value_in(self.length_units)
        y_values -= center_of_mass.y.value_in(self.length_units)
        
        sizes_of_the_markers = mass_values * 10.0
        
        plot.scatter(
            x_values, 
            y_values, 
            s=sizes_of_the_markers, 
            marker='o', 
            color="y"
        )
        
        plot.set_xlim(-4.0, 4.0)
        plot.set_ylim(-4.0, 4.0)
        

class PlotXYHistory(object):
    
    def __init__(self, indices_to_plot):
        self.length_units = nbody_system.length
        self.mass_units = nbody_system.mass
        self.indices_to_plot = indices_to_plot
        self.previous_x = [collections.deque() for x in indices_to_plot]
        self.previous_y = [collections.deque() for x in indices_to_plot]
        self.number_of_previous_points_to_show = 10
        
    def plot(self, plot, particles):
        
        if not self.indices_to_plot:
            plot.set_xlim(-1, 1)
            plot.set_ylim(-1, 1)
            return
            
        x_values = particles.x.value_in(self.length_units)
        y_values = particles.y.value_in(self.length_units)
        mass_values = particles.mass.value_in(self.mass_units)
            
        center_of_mass = particles.center_of_mass()
        
        x_values -= center_of_mass.x.value_in(self.length_units)
        y_values -= center_of_mass.y.value_in(self.length_units)
        
        for i, j in enumerate(self.indices_to_plot):
            self.previous_x[i].append(x_values[j])
            self.previous_y[i].append(y_values[j])
            
            if(len(self.previous_x[i]) > self.number_of_previous_points_to_show):
                self.previous_x[i].popleft()
                self.previous_y[i].popleft()
            
            plot.plot(self.previous_x[i],self.previous_y[i])
            
            
        
        plot.set_xlim(-2, 2)
        plot.set_ylim(-2, 2)
        del x_values
        del y_values
        
        
        
class MakeAMovie(object):
    
    
    def __init__(self, name_of_the_hdf_file):
        storage = store.StoreHDF(name_of_the_hdf_file)
        self.particles = storage.load()
        self.number_of_particles = len(self.particles)
        
    def new_list_of_indices_of_particles_to_show_in_detail(self):
        indices = []
        random_number_generator = random.Random()
        random_number_generator.seed()
        
        fraction = 40.0 / float(self.number_of_particles)
        
        for x in range(self.number_of_particles):
            if (random_number_generator.random() <= fraction):
                indices.append(x)
        return indices
    
    def new_list_of_indices_of_particles_that_stay_in_the_center(self):
        indices = set(range(self.number_of_particles))
        counter = 0
        for data in self.particles.history:
            x_values = data.x.value_in(nbody_system.length)
            y_values = data.y.value_in(nbody_system.length)
            
            center_of_mass = data.center_of_mass()
            counter += 1
            if len(indices) < 50 or counter > 200:
                break
            x_values -= center_of_mass.x.value_in(nbody_system.length)
            y_values -= center_of_mass.y.value_in(nbody_system.length)
            for i in list(indices):
                if x_values[i] < -1.0 or x_values[i] > 1.0:
                    indices.remove(i)
                elif y_values[i] < -1.0 or y_values[i] > 1.0:
                    indices.remove(i)
        print indices
        return sorted(indices)
                    
    def new_list_of_indices_of_massive_particles(self):
        mass_values = self.particles.history.next().mass
        indices_of_massive_stars = []
        for i in range(self.number_of_particles):
            if mass_values[i] > (1.0 / self.number_of_particles) | nbody_system.mass:
                indices_of_massive_stars.append(i)
        return indices_of_massive_stars
        
    def start(self):
        
        index = 0
    
        indices_of_massive_stars = []
        
        make_plot_objects = []
        make_plot_objects.append(PlotXYandMass())
        make_plot_objects.append(PlotXYandMass(self.new_list_of_indices_of_massive_particles()))
        make_plot_objects.append(PlotXYHistory(self.new_list_of_indices_of_particles_to_show_in_detail()))
        make_plot_objects.append(PlotXYHistory(self.new_list_of_indices_of_massive_particles())) #self.new_list_of_indices_of_particles_that_stay_in_the_center()))
        
        
        number_of_plots = len(list(self.particles.history))
        t0 = time.time()
        
        print "starting..."
        figure = pyplot.figure(figsize=(6,6))
        for data in self.particles.history:
        
            for i, make_a_plot in enumerate(make_plot_objects):
                plot = figure.add_subplot(2,2,i+1)
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
            "movie-cluster-10fps.avi"
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
            "movie-cluster-25fps.avi"
        ])
            
        
        #for i in range(number_of_plots):
        #    filename = "frame_"+format(i+1, "05d")+".png"
        #    os.remove(filename)
            
if __name__ == '__main__':
    use_case = MakeAMovie(sys.argv[1])
    use_case.start()
