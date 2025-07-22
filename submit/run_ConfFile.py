import sys
import os
import subprocess
import job_mananger_conf as jm
import shutil

signal_tag = 'AtoGG'

# Definitions to replace lines in executable and to generate the AOD/MiniAOD files
def replace_lines_in_file(file_path, lines_to_replace, new_lines):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    for i, line_num in enumerate(lines_to_replace):
        lines[line_num - 1] = new_lines[i] + '\n'  
    with open(file_path, 'w') as file:
        file.writelines(lines)

def generate_signal_point(
        executable,
        job_name,
        inputfile, 
        version,
        MoE,
        active_dir,
        file_path,
        condor=False
    ):

    ## Modify the executable with the specified new information
    lines_to_replace = [3, 4, 5, 6] 

    if condor:
        new_lines = ['inputfile=' + inputfile,
                    'version=' + version,
                    'MoE=' + MoE,
                    'active_directory=/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/' + active_dir] 
        execfile = file_path + executable
        replace_lines_in_file(execfile, lines_to_replace, new_lines)
        jm.submit_condor(execfile, job_name)
    else:
        new_lines = ['inputfile=' + inputfile,
                    'version=' + version,
                    'MoE=' + MoE,
                    'active_directory=/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/' + active_dir]  
        execfile = file_path + executable
        replace_lines_in_file(execfile, lines_to_replace, new_lines)
        temp = f"{file_path}{executable}"
        # temp = f"./{executable}"
        subprocess.run([temp]) 

# Setting arguments and calling event generation
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Run CLUE clustering on MiniAOD')
    parser.add_argument('--input_file', '-f', type=str, help='Specify input MiniAOD file.')
    parser.add_argument('--output_base', '-o', type=str, help='Specify the output base directory.')
    parser.add_argument('--condor', '-c', action='store_true', help='Submit jobs to condor.')
    parser.add_argument('--version', '-v', type=str, help='Label version number.')
    parser.add_argument('--MoE', '-m', type=str, help='Specify the MoE range or Ma value.')

    args = parser.parse_args()

    # Run ConfFile_cfg.py
    print("Generating CLUE Clusters")

    executable = f"run_ConfFile_ex.sh"
    outputbase = args.output_base
    job_name_temp = signal_tag + '_' + str(args.version) + 'v' + str((args.MoE)) + 'MoE' + '_CLUE'
    job_name = job_name_temp.replace(".", "p")
    file_path = f'/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/submit/'
    jobexecutable = f"{job_name}.sh"
    shutil.copy(file_path+executable, file_path+jobexecutable)
    subprocess.run(["chmod", "+x", file_path+jobexecutable], check=True)
    # file_path = f'/tmp/{USER}/DiphotonGun/'

    generate_signal_point(
        jobexecutable,
        job_name,
        str(args.input_file),
        str(args.version),
        str(args.MoE),
        str(args.output_base),
        file_path,
        condor=args.condor)