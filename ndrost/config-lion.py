#
# support/config.py.  Generated from configpy.in by configure.
#

class interpreters(object):
    python = '/opt/local/bin/python2.6' 

class compilers(object):
    cxx = 'g++-mp-4.4' 
    cc  = 'gcc-mp-4.4'
    fc = 'gfortran-mp-4.4'
    
    cxx_flags = '-m64'
    cc_flags  = '-m64'
    fc_flags = '-m64'
    
    found_fftw = 'yes'
    fftw_flags = '-I/opt/local/include  '
    fftw_libs = '-L/opt/local/lib -lfftw3 -lm  '
    
    found_gsl = 'yes'
    gsl_flags = '-I/opt/local/include  '
    gsl_libs = '-L/opt/local/lib -lgsl -lgslcblas -lm  '
    

class mpi(object):
    mpicxx = 'openmpicxx' 
    mpicc  = 'openmpicc'
    mpif95 = 'openmpif90'

class java(object):
    jni_includes = '-I/System/Library/Frameworks/JavaVM.framework/Headers'
    jdk = '/System/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Home'

class cuda(object):
    is_enabled   = 'yes'=='yes'
    compiler     = '/usr/local/cuda/bin/nvcc'
    toolkit_path = '/usr/local/cuda'
    sdk_path     = '/Developer/CUDA/C'
    
class modules(object):
    have_matplotlib = 1==1
    have_h5py       = 1==1
