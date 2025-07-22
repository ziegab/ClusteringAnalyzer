import os
import sys

import time

# preprocessing_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
# top_dir = os.path.dirname(preprocessing_dir)
# sys.path.append(top_dir)
# print(top_dir)

# import htcondor
import signal

import textwrap

# cache_dir = f"{preprocessing_dir}/cache"
# condor_cache = f"{cache_dir}/condor"

# find . -type d -exec fs setacl {} nd_campus rlidw \;

USER=os.environ['USER']
# condordir = f'/scratch365/{USER}/DiphotonGun/condor'
# condordir = f'/tmp/{USER}'
condordir = f'/afs/cern.ch/user/g/gziemyte/diphogen/cache'
condor_cache = f'{condordir}/condor'

condor_dirs = ['log', 'err', 'out']
for d in condor_dirs:
    if not os.path.isdir(f'{condordir}/{d}'):
        os.makedirs(f'{condordir}/{d}')

def submit_condor(executable, job_name, arguments=None):
    '''Submit a job to condor'''
    
    if not os.path.isfile(executable):
        # Create a dummy executable
        exec_file = f'{condor_cache}/{job_name}.sh'
        with open(exec_file, 'w') as f:
            f.write(f'#!/bin/bash\n{executable}')
        os.system(f'chmod +x {exec_file}')
        executable = exec_file

    if isinstance(arguments, list):
        arg_string = "\n            ".join(arguments)
        submit = textwrap.dedent(
            f'''\
            executable = {executable}
            log = {condordir}/log/{job_name}_$(Process).log
            MY.WantOS = "el9"
            +JobFlavour = "testmatch"

            should_transfer_files = YES
            when_to_transfer_output = ON_EXIT_OR_EVICT
            transfer_output_files = ""
            transfer_input_files = ""

            request_cpus = 1
            request_memory = 10000
            +RequestDisk = 10000000
            queue arguments from (
                {arg_string}
            )
            ''')
    else:
        submit = textwrap.dedent(
            f'''\
            executable = {executable}
            log = {condordir}/log/{job_name}.log
            MY.WantOS = "el9"
            +JobFlavour = "testmatch"

            should_transfer_files = YES
            when_to_transfer_output = ON_EXIT_OR_EVICT
            transfer_output_files = ""
            transfer_input_files = ""

            request_cpus = 1
            request_memory = 10000
            +RequestDisk = 10000000
            queue
            ''')

    submit_file = f'{condor_cache}/{job_name}.submit'
    with open(submit_file, 'w') as f:
        f.write(submit)
    
    os.system(f'condor_submit {submit_file}')
