import os

fileList = []
xrootd='root://cmsxrootd.fnal.gov/' # FNAL
#dirName = '/eos/uscms/store/group/lpcml/rchudasa/NTuples/DYToTauTau_M-50_13TeV-powheg_pythia8/DYToTauTau_ntuples-AOD/231013_070551/0000/'
dirName = '/eos/uscms/store/group/lpcml/rchudasa/NTuples/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8/WJets_ntuples_AOD/231013_055218/0000'
#dirName='/eos/uscms/store/group/lpcml/rchudasa/NTuples/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/QCD_Pt-15to7000_ntuples-AOD/231013_065542/0000/'
#dirName='/eos/uscms/store/group/lpcml/rchudasa/NTuples/TTToHadronic_TuneCP5_13TeV_powheg-pythia8/ttbar-ntuples-AOD-v2/231016_163729/0000/'
for roots,dirs,files in os.walk(dirName):
    for name in files:
	if name.endswith('root'):
		fnalDir = dirName[dirName.find('/store'):]
        	file = os.path.join(xrootd+fnalDir,name)
        	fileList.append(file)
		print(file)

#dirName2 = '/eos/uscms/store/group/lpcml/rchudasa/NTuples/DYToTauTau_M-50_13TeV-powheg_pythia8/DYToTauTau_ntuples-AOD/231013_070551/0001/'
#for roots,dirs,files in os.walk(dirName2):
#    for name in files:
#        if name.endswith('root'):
#		fnalDir2 = dirName2[dirName2.find('/store'):]
#                file = os.path.join(xrootd+fnalDir2,name)
#                fileList.append(file)

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
        #   continue 
        # Generate a unique job name for each script
        job_name = "jobConvertRootToPq_%d"%(i)
        inputFileList = ','.join(dividedFileList[i])
        #outputFile = "/eos/cms/store/group/phys_heavyions/rchudasa/e2e/ParquetFiles_correctTrackerLayerHits_SecVtxInfoAdded/DYToTauTau_M-50_13TeV-powheg_pythia8/AODJets-inference/DYToTauTau_M-50_13TeV-powheg_pythia8_%d.parquet"%(i)
        outputFile = "/eos/cms/store/group/phys_heavyions/rchudasa/e2e/ParquetFiles_correctTrackerLayerHits_SecVtxInfoAdded/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8/AODJets-Inference/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8_%d.parquet"%(i)
        #outputFile = "/eos/cms/store/group/phys_heavyions/rchudasa/e2e/ParquetFiles/DYToTauTau_M-50_13TeV-powheg_pythia8/AODJets/DYToTauTau_M-50_13TeV-powheg_pythia8_%d.parquet"%(i)
        #outputFile = "/eos/cms/store/group/phys_heavyions/rchudasa/e2e/ParquetFiles/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8/AODJets/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8_%d.parquet"%(i)
        #outputFile = "/eos/cms/store/group/phys_heavyions/rchudasa/e2e/ParquetFiles_correctTrackerLayerHits_SecVtxInfoAdded/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/AODJets-Inference/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8_%d.parquet"%(i)
        #outputFile = "/eos/cms/store/group/phys_heavyions/rchudasa/e2e/ParquetFiles_correctTrackerLayerHits_SecVtxInfoAdded/TTToHadronic_TuneCP5_13TeV_powheg-pythia8/AODJets-Inference/TTToHadronic_TuneCP5_13TeV_powheg-pythia8_%d.parquet"%(i)
        #print(type(outputFile))
        # Define the contents of the SLURM bash script
        script_contents = """#!/bin/bash
export X509_USER_PROXY=$1
voms-proxy-info -all
voms-proxy-info -all -file $1
source /cvmfs/sft.cern.ch/lcg/views/LCG_97a/x86_64-centos7-gcc8-opt/setup.sh
source /afs/cern.ch/work/r/rchudasa/private/venv/bin/activate
cd /afs/cern.ch/work/r/rchudasa/private/TauClassification/miniAOD_checks/MLAnalyzer/convertRootFiles/
echo $PWD
python /afs/cern.ch/work/r/rchudasa/private/TauClassification/miniAOD_checks/MLAnalyzer/convertRootFiles/convert_root2pq_tau_jet.py -i %s -o %s
"""%(inputFileList,outputFile)
	
	condor_contents = """Proxy_path            = /afs/cern.ch/work/r/rchudasa/private/x509up_u43677 
arguments             = $(Proxy_path) $(ProcId)
executable            = %s.sh
Universe              = vanilla
GetEnv                = True
output                = output/$(ClusterId).$(ProcId).out
error                 = error/$(ClusterId).$(ProcId).err
log                   = log/$(ClusterId).log
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
#output_directory = "condorScriptsTTbar_Inference"
#output_directory = "condorScriptsQCD_Inference"
output_directory = "condorScriptsWJets_Inference"
#output_directory = "condorScriptsDYTauTau_Inference"
#output_directory = "condorScriptsDYTauTau"
#output_directory = "condorScriptsDYEE"
#output_directory = "condorScriptsQCDEMEnriched"

generate_condor_scripts(len(dividedFileList), output_directory)
