#!/usr/bin/env python

import os,sys
import argparse
import commands
import math
import random
import subprocess

random.seed()

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("-m", "--model"     , dest="model"       , help="set detector model (0-3)"      , default=2, type=int)
parser.add_argument("-v", "--version"   , dest="version"     , help="set detector version (1-7)"    , default=1, type=int)
parser.add_argument("-b", "--Bfield"    , dest="Bfield"      , help="B field value in Tesla"        , default=0, type=float)
parser.add_argument("-f", "--lhefile"   , dest="lhefile"     , help="name of file in directory"     , required=True)
parser.add_argument("-j", "--numjobs"   , dest="numjobs"     , help="number of jobs to run"         , default=1, type=int)
parser.add_argument("-o", "--dataOutDir", dest="dataOutDir"  , help="directory to output root files", required=True)
parser.add_argument("-S", "--no-submit" , action="store_true", dest="nosubmit" , help="Do not submit batch job.")
arg = parser.parse_args()

dataOutDir = arg.dataOutDir
# Check for trailing slash on lhedir and outdir and delete
if arg.dataOutDir.split("/")[-1] == "": dataOutDir = arg.dataOutDir[:-1]

filename = arg.lhefile.split("/")[-1]
nevts = int(filename.split("_")[0])
outFilename = str(filename.split(".lhe")[0])

# Create subdirectory based on particle type
particle = filename.split("_")[-1].split(".lhe")[0]

# Check that the lhe and output directories exist
if not os.path.exists(dataOutDir):
    print "Provided output directory \"%s\" does not exist!"%(dataOutDir)
    quit()

# Check for temp directory and create one if not
if not os.path.exists("./temp"): os.mkdir("temp")

bval="BOFF"
if arg.Bfield>0 : bval="BON" 

outDir = os.getcwd()
outTag="version%d_model%d_%s"%(arg.version,arg.model,outFilename)

# Wrapper
scriptFile = open("%s/runJob.sh"%(outDir), "w")
scriptFile.write("#!/bin/bash\n")
scriptFile.write("job=$1\n")
scriptFile.write("name=$2\n")
scriptFile.write("source /data/cmszfs1/sw/HGCAL_SIM_A/setup.sh\n")
scriptFile.write("cd temp\n")
scriptFile.write("mkdir ${name}_${job}\n")
scriptFile.write("cd ${name}_${job}\n")
scriptFile.write("cp $PWD/../../g4steer_${name}_${job}.mac g4steer_${name}_${job}.mac\n")
scriptFile.write("rm $PWD/../../g4steer_${name}_${job}.mac\n")
scriptFile.write("./../../PFCalEE g4steer_${name}_${job}.mac %d %d 1 | tee g4.log\n"%(arg.version,arg.model))
scriptFile.write("mv PFcal.root %s/HGcal_%s_${job}.root\n"%(dataOutDir,outTag))
scriptFile.write("localdir=`pwd`\n")
scriptFile.write("echo \"--Local directory is \" $localdir >> g4.log\n")
scriptFile.write("ls * >> g4.log\n")
scriptFile.write("cd ..\n")
scriptFile.write("rm -r ${name}_${job}\n")
scriptFile.write("cd ../..\n")
scriptFile.write("echo \"All done\"\n")
scriptFile.close()

# Submit to condor
condorSubmit = open("%s/condorSubmit"%(outDir), "w")
condorSubmit.write("Executable          =  %s\n" % scriptFile.name)
condorSubmit.write("Universe            =  vanilla\n")
condorSubmit.write("Requirements        =  Arch==\"X86_64\"  &&  (Machine  !=  \"zebra01.spa.umn.edu\")  &&  (Machine  !=  \"zebra02.spa.umn.edu\")  &&  (Machine  !=  \"zebra03.spa.umn.edu\")  &&  (Machine  !=  \"zebra04.spa.umn.edu\")\n")
condorSubmit.write("+CondorGroup        =  \"cmsfarm\"\n")
condorSubmit.write("getenv              =  True\n")
condorSubmit.write("Request_Memory      =  4 Gb\n")
condorSubmit.write("Log         =  %s.log\n" % outDir)


for job in xrange(0,arg.numjobs):

    nevts = int(filename.split("_")[0])
    outFilename = str(filename.split(".lhe")[0])
    
    # Finish writing condor file
    condorSubmit.write("Arguments       = %d %s\n"%(job,outTag))
    condorSubmit.write("Queue\n")

    # Write GEANT4 macro
    g4Macro = open("%s/g4steer_%s_%d.mac"%(outDir,outTag,job), "w")
    g4Macro.write("/control/verbose 0\n")
    g4Macro.write("/run/verbose 0\n")
    g4Macro.write("/event/verbose 0\n")
    g4Macro.write("/tracking/verbose 0\n")
    g4Macro.write("/N03/det/setField %1.1f T\n"%(arg.Bfield))
    g4Macro.write("/random/setSeeds %d %d\n"%(random.uniform(0,100000000),random.uniform(0,100000000)))
    g4Macro.write("/filemode/inputFilename %s\n"%(arg.lhefile))
    g4Macro.write("/run/initialize\n")
    g4Macro.write("/run/beamOn %d\n"%(nevts))
    g4Macro.close()

condorSubmit.close()

os.system("chmod u+rwx %s/runJob.sh"%(outDir))

command = "condor_submit " + condorSubmit.name + "\n"
subprocess.call(command.split())

