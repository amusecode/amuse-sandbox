# -*- coding: utf-8 -*-
"""
Utility routines for doing computational experiments.
Jürgen Jänes <jurgen.janes@gmail.com>

Warning: This package is mscware and needs significant refactoring and documentation
to become effortlessly usable.
"""

import commands
import datetime
import math
import os
import re
import pickle
import shutil
import subprocess
import sys

import __main__

import numpy as np
#import scipy as sp

# calculates the mean and uncertainty from a list of measurements
# assumes standard distribution and 95% confidence intervals
def sprintu(values, format):

    avg = sp.mean(values)
    err = 1.96 * sp.std(values) / math.sqrt(len(values)) 
    fmt = "%s±%s" % (format, format)
    return fmt % (avg, err)

# runs an external command and returns the output or fails gracefully
def run_or_fail(s):

    print "running: %s" % (s, )
    r = commands.getstatusoutput( s )
    if (r[0] != 0):
        print "error executing: %s" % (s, )
        print "output: %s" % (r[1], )
        exit(0)
    return(r)

# calculate row means
def means(a):

    m = []
    for i in range(len(a)):
        m.append(sp.mean(a[i]))

    return m

# calculate row errors
def errors(a):

    m = []
    for i in range(len(a)):
        m.append(1.96 * sp.std(a[i]) / math.sqrt(len(a[i])))

    return m

# create a directory if it does not exist
def mkresdir(exp_name):
    # TODO add microsecond-resolutions?
    timestamp = datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    dirname = "%s-%s" % (exp_name, timestamp)
    if not os.path.isdir(dirname):
        os.mkdir( dirname )
    return (dirname)
    #print mkresdir("compare_models")

def mkresfile(exp_name, file_name):
    return os.path.join(exp_name, file_name)

def save_variable(val, exp_name, file_name):
    #cr.save_variable((self.x1, self.y1, self.x2, self.y2), dout, "%d-%d-coordinates" % (iSlice, iROI), v)
    outf = open(mkresfile(exp_name, '%s.pkl' % (file_name)), 'wb')
    pickle.dump(val, outf)
    outf.close()
    print "save_variable: wrote to %s" % (outf,)

def load_variable(exp_name, file_name):
    #(self.x1, self.y1, self.x2, self.y2) = cr.load_variable(din, "%d-%d-coordinates" % (iSlice, iROI), v)
    inpf = open(mkresfile(exp_name, '%s.pkl' % (file_name)), 'rb')
    val = pickle.load(inpf)
    inpf.close()
    print "load_variable: read from %s" % (inpf,)
    return val

class ResultsWriter:
    def __init__(self, exp_name=None):
        # use the name of the main script file as the name of the output directory
        if (exp_name is None):
            exp_name = os.path.splitext(__main__.__file__)[0]

        self.exp_name = exp_name
        self.exp_fname = mkresdir(exp_name)

        # ugly hack to duplicate screen output to a file
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
        tee = subprocess.Popen(["tee", "%s/script-output.txt" % (self.exp_fname)], stdin=subprocess.PIPE)
        os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
        os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

        # make a copy of the s    
        print "ic: making copy of the driver script\n\tfrom: %s\n\tto: %s" % (__main__.__file__, self.exp_fname)
        shutil.copy(__main__.__file__, self.exp_fname)

        self.computation_start = datetime.datetime.now()
        print "computation start: %s" % (self.computation_start, )

    def __del__(self):
        # ugly hack (otherwise ResultsReader will fail)
        try:
            computation_stop = datetime.datetime.now()
            print "computation start: %s" % (self.computation_start, )
            print "computation stop: %s" % (computation_stop, )
            print "total running time: %s" % (computation_stop - self.computation_start, )
        except AttributeError:
            pass

    def fn(self, file_name):
        return mkresfile(self.exp_fname, file_name)
    
    def save_variable(self, val, file_name):
        save_variable(val, self.exp_fname, file_name)

    def cdto_expdir(self, subDir = None):
        self.prevdir = os.getcwd()
        newdir = self.exp_fname if (subDir is None) else os.path.join(self.prevdir, self.exp_fname, subDir)
        if not os.path.isdir(newdir):
            os.mkdir(newdir)
        os.chdir(newdir)
        print "ic: now in directory: %s" % (newdir,)

    def cdto_prevdir(self):
        os.chdir(self.prevdir)
        print "ic: now in directory: %s" % (self.prevdir,)

class ResultsReader(ResultsWriter):
    def __init__(self, exp_fname=None):
        self.exp_name = os.path.basename( os.path.splitext(__main__.__file__)[0] )
        if not (exp_fname is None):
            self.exp_fname = exp_fname
        else:
            results = filter(os.path.isdir, os.listdir("./"))
            results.sort()
            p = re.compile("^%s" % self.exp_name)
            results = filter(p.match, results)
            self.exp_fname = results[-1]
            print "setting exp_fname to: %s" % (self.exp_fname)
        print "ic: loading results from: %s" % (self.exp_fname)

    def load_variable(self, file_name):
        return load_variable(self.exp_fname, file_name)
    
if __name__ == '__main__':
    exp = mkresdir("initialconditions_test")
    a = [1, 2, 3]
    b = "Aloha!"
    
    save_variable((a, b), exp, "a_and_b")
    (a_new, b_new) = load_variable(exp, "a_and_b")
    
    print a_new
    print b_new
