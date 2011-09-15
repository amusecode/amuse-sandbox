import subprocess
import os.path
import imp
import sys
import hashlib
import matplotlib
import matplotlib.pyplot as plt
import shutil
import cStringIO
import numpy.random
import time
import json
import logging

from matplotlib import _pylab_helpers
from amuse.community import get_amuse_root_dir

             
def clear_state():
    plt.close('all')
    matplotlib.rc_file_defaults()
    
def run_code(plot_path):
    """
    Import a Python module from a path, and run the function given by
    name, if function_name is not None.
    """
    pwd = os.getcwd()
    path, fname = os.path.split(plot_path)
    sys.path.insert(0, os.path.abspath(path))
    stdout = sys.stdout
    #sys.stdout = cStringIO.StringIO()
    os.chdir(path)
    fd = None
    try:
        fd = open(fname)
        print fname
        numpy.random.seed(1234)
        module = imp.load_module(
            "__plot__", fd, fname, ('py', 'r', imp.PY_SOURCE))
    finally:
        del sys.path[0]
        os.chdir(pwd)
        sys.stdout = stdout
        if fd is not None:
            fd.close()

def run_savefig(plot_path, basename, tmpdir, formats):
    """
    Once a plot script has been imported, this function runs savefig
    on all of the figures in all of the desired formats.
    """
    fig_managers = _pylab_helpers.Gcf.get_all_fig_managers()
    
    outpaths = []
    
    for i, figman in enumerate(fig_managers):
        for (format, dpi) in formats:
            if len(fig_managers) == 1:
                outname = basename
            else:
                outname = "%s_%02d" % (basename, i)
            outname = outname + "." + format
            outpath = os.path.join(tmpdir, outname)
            try:
                figman.canvas.figure.savefig(outpath, dpi=dpi)
            except:
                s = cbook.exception_to_str("Exception saving plot %s" % plot_path)
                warnings.warn(s, PlotWarning)
                return 0
            outpaths.append(outpath)
            
    return outpaths
        
            
def render_figures(plot_path, tmpdir, formats):
    """
    Run a pyplot script and save the low and high res PNGs and a PDF
    in outdir.
    """
    plot_path = str(plot_path)  # todo, why is unicode breaking self
    basedir, fname = os.path.split(plot_path)
    basename, ext = os.path.splitext(fname)

    clear_state()
    
    run_code(plot_path)
    
    names = run_savefig(plot_path, basename, tmpdir, formats)

    if '__plot__' in sys.modules:
        del sys.modules['__plot__']

    return names





class CheckExamples(object):

    def __init__(self):
        simple_examples_dir = os.path.join(get_amuse_root_dir(), 'examples', 'simple')
        self.simple_scripts = [os.path.join(simple_examples_dir, one_file) for one_file in os.listdir(simple_examples_dir) if one_file[-3:] == '.py']
        
        
    def run_example(self, path):
        logging.disable(logging.INFO)
        t0 = time.time()
        names = render_figures(path, 'tmp', [('png', 80),])
        t1 = time.time()
        return names, t1 - t0
    
    
    def make_references(self):
        for x in self.simple_scripts:
            scriptname = os.path.basename(x)
            if scriptname.startswith('_'):
                continue
            figures, time = self.run_example(x)
            digests = [self.get_digest(path) for path in figures]
            yield [os.path.basename(x), time, digests]
                
    
    def get_digest(self, path):
        with open(path, 'rb') as figurefile:
            allbytes = figurefile.read()
        return hashlib.sha256(allbytes).hexdigest()
   
   
if __name__ == '__main__':
    uc = CheckExamples()
    all = []
    for x in uc.make_references():
        print x
        all.append(x)
    print json.dumps(all, indent = 4)
    
