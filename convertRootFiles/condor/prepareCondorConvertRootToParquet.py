import os

fileList = []
for roots,dirs,files in os.walk('/eos/uscms/store/group/lpcml/rchudasa/NTuples/DYToTauTau_M-50_13TeV-powheg_pythia8/DYToTauTau_ntuples_miniAOD/231005_120502/0000/'):
    for name in files:
	if name.endswith('root'):
        	file = os.path.join(roots,name)
        	fileList.append(file)

for roots,dirs,files in os.walk('/eos/uscms/store/group/lpcml/rchudasa/NTuples/DYToTauTau_M-50_13TeV-powheg_pythia8/DYToTauTau_ntuples_miniAOD/231005_120502/0001/'):
    for name in files:
        if name.endswith('root'):
                file = os.path.join(roots,name)
                fileList.append(file)

def divide_list(big_list, chunk_size):
    divided_lists = []
    for i in range(0, len(big_list), chunk_size):
        chunk = big_list[i:i + chunk_size]
        divided_lists.append(chunk)
    return divided_lists

dividedFileList = divide_list(fileList,15)

def generate_condor_scripts(numScripts, output_directory):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
	os.makedirs(output_directory)

    for i in range(numScripts):
        #if i < 2:
        #    continue
        # Generate a unique job name for each script
        job_name = "jobMergeParquet_%d"%(i)
        inputFileList = ','.join(dividedFileList[i])
        outputFile = "/eos/uscms/store/group/lpcml/rchudasa/ParquetFiles/DYToTauTau_M-50_13TeV-powheg_pythia8/miniAODJets/DYToTauTau_M-50_13TeV-powheg_pythia8_%d.parquet"%(i)
        #print(type(outputFile))
        # Define the contents of the SLURM bash script
        script_contents = """#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_97a/x86_64-centos7-gcc8-opt/setup.sh
source /uscms/home/rchudasa/nobackup/venv/bin/activate
cd /uscms/home/rchudasa/nobackup/miniAOD_checks/CMSSW_10_6_25/src/MLAnalyzer/convertRootFiles/condor/
echo $PWD
python /uscms/home/rchudasa/nobackup/miniAOD_checks/CMSSW_10_6_25/src/MLAnalyzer/convertRootFiles/convert_root2pq_tau_jet.py -i %s -o %s
"""%(inputFileList,outputFile)
	
	condor_contents = """Universe              = vanilla
executable            = %s.sh
arguments             = $(ProcId)
GetEnv                = True
output                = output/$(ClusterId).$(ProcId).out
error                 = error/$(ClusterId).$(ProcId).err
requirements          = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "tomorrow"

queue 1
"""%(job_name)


        # Write the script contents to a file
        script_path = os.path.join(output_directory, "%s.sh"%(job_name))
        with open(script_path, "w") as script_file:
            script_file.write(script_contents)

        print("Generated condor script: %s"%(script_path))

        # Write the condor script 
        condor_path = os.path.join(output_directory, "%s.sub"%(job_name))
        with open(condor_path, "w") as condor_file:
            condor_file.write(condor_contents)

        print("Generated condor sub script: %s"%(condor_path))

# Example usage
#output_directory = "condorScriptsQCD"
#output_directory = "condorScriptsWJets"
output_directory = "condorScriptsDYTauTau"
#output_directory = "condorScriptsDYEE"
#output_directory = "condorScriptsQCDEMEnriched"

generate_condor_scripts(len(dividedFileList), output_directory)
