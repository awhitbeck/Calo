#!/usr/bin/env python

import os,sys
import argparse
import commands
import random
import subprocess
from time import strftime

myTime = strftime("%H%M%S")
random.seed()

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("-m", "--model"     , dest="model"       , help="set detector model (0-3)"      , default=2, type=int)
parser.add_argument("-v", "--version"   , dest="version"     , help="set detector version (1-7)"    , default=1, type=int)
parser.add_argument("-f", "--lhefile"   , dest="lhefile"     , help="full path to input file"       , required=True)
parser.add_argument("-j", "--numjobs"   , dest="numjobs"     , help="number of jobs to run"         , default=1, type=int)
parser.add_argument("-o", "--outputDir" , dest="outputDir"   , help="directory to output root files", required=True)
arg = parser.parse_args()

outputDir = arg.outputDir

# Check that the output directory exist
if not os.path.exists(outputDir):
    print "Provided output directory \"%s\" does not exist!"%(outputDir)
    print "Exiting..."
    quit()

# Check that the input .lhe file exists
if not os.path.exists(arg.lhefile):
    print "Provided input .lhe file \"%s\" does not exist!"%(arg.lhefile)
    print "Exiting..."
    quit()

# Check for trailing slash on outdir and delete
if arg.outputDir.split("/")[-1] == "": outputDir = arg.outputDir[:-1]

filename    = arg.lhefile.split("/")[-1]
nevts       = int(filename.split("_")[0])
outFilename = str(filename.split(".lhe")[0])

# Check for existence of temp directory and create one if none is found
if not os.path.exists("./temp"): os.mkdir("temp")

outDir = os.getcwd()
outTag = "version%d_model%d_%s"%(arg.version,arg.model,outFilename)

# Write .sh script to be run by Condor 
scriptFile = open("%s/runJob_%s.sh"%(outDir,myTime), "w")
scriptFile.write("#!/bin/bash\n")
scriptFile.write("job=$1\n")
scriptFile.write("name=$2\n")
scriptFile.write("source /data/cmszfs1/sw/HGCAL_SIM_A/setup.sh\n")
scriptFile.write("cd temp\n")
scriptFile.write("mkdir ${name}_${job}\n")
scriptFile.write("cd ${name}_${job}\n")
scriptFile.write("cp $PWD/../../g4steer_${name}_${job}.mac g4steer.mac\n")
scriptFile.write("cp $PWD/../../b18d36.dat b18d36.dat\n")
scriptFile.write("rm $PWD/../../g4steer_${name}_${job}.mac\n")
scriptFile.write("./../../PFCalEE g4steer.mac %d %d 0\n"%(arg.version,arg.model))
scriptFile.write("mv PFcal.root %s/HGcal_%s_${job}.root\n"%(outputDir,outTag))
scriptFile.write("cd ..\n")
scriptFile.write("rm -r ${name}_${job}\n")
scriptFile.write("cd ../..\n")
scriptFile.write("echo \"All done\"\n")
scriptFile.close()

# Write Condor submit file 
condorSubmit = open("%s/condorSubmit_%s"%(outDir,myTime), "w")
condorSubmit.write("Executable          =  %s\n" % scriptFile.name)
condorSubmit.write("Universe            =  vanilla\n")
condorSubmit.write("Requirements        =  Arch==\"X86_64\"  &&  (Machine  !=  \"zebra01.spa.umn.edu\")  &&  (Machine  !=  \"zebra02.spa.umn.edu\")  &&  (Machine  !=  \"zebra03.spa.umn.edu\")  &&  (Machine  !=  \"zebra04.spa.umn.edu\")\n")
condorSubmit.write("+CondorGroup        =  \"cmsfarm\"\n")
condorSubmit.write("getenv              =  True\n")
condorSubmit.write("Request_Memory      =  4 Gb\n")
condorSubmit.write("Log                 =  %s.log\n" % outDir)

for job in xrange(10000,10000+arg.numjobs):

    # Append jobs to Condor submit file 
    condorSubmit.write("Arguments       = %d %s\n"%(job,outTag))
    condorSubmit.write("Queue\n")

    # Write GEANT4 macro for each job
    g4Macro = open("%s/g4steer_%s_%d.mac"%(outDir,outTag,job), "w")
    g4Macro.write("/control/verbose 0\n")
    g4Macro.write("/run/verbose 0\n")
    g4Macro.write("/event/verbose 0\n")
    g4Macro.write("/tracking/verbose 0\n")
    g4Macro.write("/random/setSeeds %d %d\n"%(random.uniform(0,100000000),random.uniform(0,100000000)))
    g4Macro.write("/filemode/inputFilename %s\n"%(arg.lhefile))
    g4Macro.write("/run/initialize\n")
    g4Macro.write("/run/beamOn %d\n"%(nevts))
    g4Macro.close()

condorSubmit.close()

os.system("chmod u+rwx %s/runJob_%s.sh"%(outDir,myTime))

command = "condor_submit " + condorSubmit.name + "\n"
subprocess.call(command.split())
