import subprocess
import os.path
import imp
import sys
import hashlib
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shutil
import cStringIO
import numpy.random
import time
import json
import logging
import datetime
from optparse import OptionParser

from matplotlib import _pylab_helpers
from amuse.community import get_amuse_root_dir
from amuse.support.thirdparty.texttable import Texttable
from amuse.support import literature

             
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
    sys.argv = sys.argv[:1]
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
            
            figman.canvas.figure.savefig(outpath, dpi=dpi)
            
            outpaths.append(outpath)
            
    return outpaths
        
            
def render_figures(plot_path, tmpdir, formats):
    """
    Run a pyplot script and save the low and high res PNGs and a PDF
    in outdir.
    """
    plot_path = str(plot_path)
    basedir, fname = os.path.split(plot_path)
    basename, ext = os.path.splitext(fname)

    names = run_savefig(plot_path, basename, tmpdir, formats)
    
    if '__plot__' in sys.modules:
        del sys.modules['__plot__']

    return names





class CheckExamples(object):

    def __init__(self, references_directory, current_directory):
        self.references_directory = references_directory
        self.current_directory = current_directory
        simple_examples_dir = os.path.join(get_amuse_root_dir(), 'examples', 'simple')
        self.simple_scripts = [os.path.join(simple_examples_dir, one_file) for one_file in os.listdir(simple_examples_dir) if one_file[-3:] == '.py']
        
    def run_example(self, path, directory):
        logging.disable(logging.INFO)
        t0 = time.time()
        error = ''
        clear_state()
        try:
            run_code(path)
        except Exception as ex:
            error = 'error in running code: ' + str(ex)
        try:
            names = render_figures(path, directory, [('png', 80),])
        except Exception as ex:
            error = 'error in making plot: ' + ex
            
        t1 = time.time()
        return names, t1 - t0, error
    
    
    def make_references(self):
        for x in self.simple_scripts:
            scriptname = os.path.basename(x)
            if scriptname.startswith('_'):
                continue
            figures, time, error = self.run_example(x, self.references_directory)
            digests = [self.get_digest(path) for path in figures]
            yield [os.path.basename(x), time, figures, digests, error]
               
    
    def make_current(self):
        for x in self.simple_scripts:
            scriptname = os.path.basename(x)
            if scriptname.startswith('_'):
                continue
            figures, time, error = self.run_example(x, self.current_directory)
            digests = [self.get_digest(path) for path in figures]
            yield [os.path.basename(x), time, figures, digests, error]
                 
    
    def get_digest(self, path):
        with open(path, 'rb') as figurefile:
            allbytes = figurefile.read()
        return hashlib.sha256(allbytes).hexdigest()
   
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "--checkdir", 
        default = "check",
        dest="references_directory",
        help="directory of the reference pictures",
        type="string"
    )
    result.add_option(
        "--currentdir", 
        default = "current",
        dest="current_directory",
        help="directory of the pictures to generate in this run",
        type="string"
    )
    result.add_option(
        "-r", "--generate-references",
        default = False,
        dest="must_make_references",
        help="generate the reference pictures",
        action="store_true",
    )
    return result

def get_reference(name, references):
    
    for refname, dt, figures, digests, error in references:
        if refname == name:
            return refname, dt, figures, digests, error
    return None, None, None, None


def print_table(rows):
    rows = [ [name, print_time(t1), print_time(t2), comment] for name, t1, t2, comment in rows]
    table = Texttable()
    table.set_cols_align(["l", "r", "r", "l"])
    table.set_cols_dtype(['t', 'f', 'f', 't'])  
    rows.insert(0,['Name', 'reference\ntime', 'current\ntime', 'comment'])
    table.add_rows(rows)
    print table.draw() + "\n"

def print_time(seconds):
    return  str(datetime.timedelta(seconds=round(seconds,1)))
    
def main(references_directory = 'check', current_directory='currrent', must_make_references = False):
    if not os.path.exists(references_directory):
        os.mkdir(references_directory)
    if not os.path.exists(current_directory):
        os.mkdir(current_directory)
    
    if not os.path.exists('references.json'):
        must_make_references = True
        
    uc = CheckExamples(references_directory, current_directory)
    references = []
    if must_make_references:
        for x in uc.make_references():
            references.append(x)
        with open('references.json', 'w') as stream:
            json.dump(references,stream, indent = 4)
    else:
        with open('references.json', 'r') as stream:
            references = json.load(stream)
            
    current = []
    for x in uc.make_current():
        current.append(x)
    rows = []
    has_a_difference = False
    for name, dt, figures, digests, error in current:
        refname, refdt, reffigures, refdigests, referror = get_reference(name, references)
        if refname is None:
            has_a_difference = True
            rows.append([ name, 0.0, dt, 'example not in reference']) 
        elif referror:
            has_a_difference = True
            if error == referror:
                rows.append([ name, refdt, dt,  'error in reference+check: '+referror])
            else:
                rows.append([ name, refdt, dt,  'error in reference: '+referror])
        elif error:
            has_a_difference = True
            rows.append([ name, refdt, dt,  'error in current: '+error]) 
        else:
            lookup = {}
            for reffigure, refdigest in zip(reffigures, refdigests):
                reffigure = os.path.basename(reffigure)
                lookup[reffigure] = refdigest
            print lookup
            compare = []
            for figure, digest in zip(figures, digests):
                figure = os.path.basename(figure)
                if not figure in lookup:
                    compare.append(figure  + ' not in reference')
                else:
                    refdigest = lookup[figure]
                    if not refdigest == digest:
                        compare.append(os.path.basename(figure) + ' differes from reference')
            if len(compare) > 0:
                has_a_difference = True
            else:
                compare.append('equal')
                
            rows.append([ name, refdt, dt, '\n'.join(compare)]) 
    print_table(rows)
    if has_a_difference:
        return 1
    else:
        return 0
    
   
if __name__ == '__main__':
    literature.TrackLiteratureReferences.suppress_output()
    options, arguments = new_option_parser().parse_args()
    sys.exit(main(**options.__dict__))
