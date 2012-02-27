import os
import os.path
import sys
import cProfile
import subprocess
import shutil
import time

"""
amuse.sh  -m cProfile -o run1.prof particles_and_gas_in_cluster.py   \
    --nstar=1000 \
    --end-time=1 \
    --gas-fraction=0.90 \
    --seed=1234\
    --star_smoothing_fraction=0.05 \
    --gas_smoothing_fraction=0.05 \
    --interaction-timestep=0.01 \
    --ntimesteps=10 \
    --noplot

"""
def run_bridge():
    
    script_name = '../../../examples/validation/particles_and_gas_in_cluster.py'
    common_arguments = [
        sys.executable,
        '-m', 'cProfile',
        '-o', 'run.prof',
        script_name,
        '--end-time=1',
        '--seed=1234',
        #'--interaction-timestep=0.02',
        '--noplot'
    ]
    
    stats_arguments = [
        sys.executable,
        '../show_validation_stats.py',
        'run.prof'
    ]
    environment = {}
    environment.update(os.environ)
    environment['PYTHONPATH'] = '../../../src'
    
    for n in [500, 1000, 10000]:
        for starcode in ['fi', 'octgrav', 'phigrape' ]:
            for gascode in ['fi', 'field']:
		print "start run , starcode:", starcode, ", gascode:", gascode, " n:", n
                runname = 'run2_star_{0}__gas_{1}__n_{2}'.format(starcode, gascode, n)
                
                if os.path.exists(runname):
                    is_backuped = False
                    for x in range(10):
                        backupname = runname + '.{0}'.format(x)
                        if os.path.exists(backupname):
                            continue
                        os.rename(runname, backupname)
                        is_backuped = True
                        break
                    if not is_backuped:
                        shutil.rmtree(backupname+'.0')
                        os.rename(runname, backupname)
                os.mkdir(runname)
                arguments = list(common_arguments)
                arguments.append('--gas-code={0}'.format(gascode))
                arguments.append('--star-code={0}'.format(starcode))
                arguments.append('--nstar={0}'.format(n))
                
                with open(os.path.join(runname, 'run.out'), 'wb') as outfile:
                    process = subprocess.Popen(
                        arguments,
                        cwd = runname, 
                        stdout = outfile,
                        stderr = subprocess.STDOUT,
                        env = environment
                    )
                
                    process.wait()
                
                
                time.sleep(4.0)
                
                with open(os.path.join(runname, 'run.stat'), 'wb') as outfile:
                    process = subprocess.Popen(
                        stats_arguments,
                        cwd = runname, 
                        stdout = outfile,
                        stderr = subprocess.STDOUT,
                        env = environment
                    )
                    process.wait()
		print "end run , starcode:", starcode, ", gascode:", gascode, " n:", n
                    

def run_sph():
    
    script_name = '../../../examples/validation/particles_and_gas_in_cluster.py'
    common_arguments = [
        sys.executable,
        '-m', 'cProfile',
        '-o', 'run.prof',
        script_name,
        '--end-time=1',
        '--seed=1234',
        '--allinoone',
        '--noplot'
    ]
    
    stats_arguments = [
        sys.executable,
        '../show_validation_stats.py',
        'run.prof'
    ]
    environment = {}
    environment.update(os.environ)
    environment['PYTHONPATH'] = '../../../src'
    
    for n in [500, 1000, 10000]:
        for sphcode in ['fi' ]:
            runname = 'run2_sph_{0}__n_{1}'.format(sphcode, n)
            
            if os.path.exists(runname):
                is_backuped = False
                for x in range(10):
                    backupname = runname + '.{0}'.format(x)
                    if os.path.exists(backupname):
                        continue
                    os.rename(runname, backupname)
                    is_backuped = True
                    break
                if not is_backuped:
                    shutil.rmtree(backupname+'.0')
                    os.rename(runname, backupname)
            os.mkdir(runname)
            arguments = list(common_arguments)
            arguments.append('--sph-code={0}'.format(sphcode))
            arguments.append('--nstar={0}'.format(n))
            
            with open(os.path.join(runname, 'run.out'), 'wb') as outfile:
                process = subprocess.Popen(
                    arguments,
                    cwd = runname, 
                    stdout = outfile,
                    stderr = subprocess.STDOUT,
                    env = environment
                )
            
                process.wait()
            
            
            time.sleep(4.0)
            
            with open(os.path.join(runname, 'run.stat'), 'wb') as outfile:
                process = subprocess.Popen(
                    stats_arguments,
                    cwd = runname, 
                    stdout = outfile,
                    stderr = subprocess.STDOUT,
                    env = environment
                )
                process.wait()

if __name__ == '__main__':
    run_sph()
    run_bridge()
                
