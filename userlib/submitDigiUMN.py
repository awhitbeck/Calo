#! /usr/bin/env python

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
parser.add_argument("-m", "--model"     ,    dest="model"      , help="detector model"                    , default=0, type=int)
parser.add_argument("-d", "--inputDir" ,    dest="inputDir"    , help="directory containing files"        , required=True)
parser.add_argument("-o", "--outputDir",    dest="outputDir" , help="directory to output digiroot files", required=True)
arg = parser.parse_args()

inputDir = arg.inputDir
outputDir = arg.outputDir

# Check that the root and output directories exist
if not os.path.exists(inputDir):
    print "Provided root directory \"%s\" does not exist!"%(inputDir)
    print "Exiting..."
    quit()

if not os.path.exists(outputDir):
    print "Provided output directory \"%s\" does not exist!"%(outputDir)
    print "Exiting..."
    quit()

# Check for trailing slash on inputDir and outdir and delete
if arg.inputDir.split("/")[-1] == "": inputDir = arg.inputDir[:-1]
if arg.outputDir.split("/")[-1] == "": outputDir = arg.outputDir[:-1]

# Check for temp directory and create one if none exists 
if not os.path.exists("./temp"): os.mkdir("temp")

# Digitization parameters 
# Format for granularity, noise, threshold: beginLayer-endLayer:value,beginLayer2-endLayer2:value2...
granularity = "0-55:4"
noise       = "0-55:0.15" # In MIPs
threshold   = "0-55:4"    # In ADC
interCalib  = 3           # In %

outDir = os.getcwd()

# Write .sh file to be submitted to Condor
scriptFile = open("%s/runDigiJob.sh"%(outDir), "w")
scriptFile.write("#!/bin/bash\n")
scriptFile.write("infile=$1\n")
scriptFile.write("outfilepath=$2\n")
scriptFile.write("outfilename=$3\n")
scriptFile.write("source /data/cmszfs1/sw/HGCAL_SIM_A/setup.sh\n")
scriptFile.write("localdir=`pwd`\n")
scriptFile.write("./bin/digitizer -N 0 -O ${outfilepath} -F ${outfilename} -I ${infile} $PWD -G %s -S %s -T %s -C %d\n"%(granularity,noise,threshold,interCalib))
scriptFile.write("echo \"All done\"\n")
scriptFile.close()

# Write Condor submit file
condorSubmit = open("%s/condorSubmit"%(outDir), "w")
condorSubmit.write("Executable          =  %s\n"%(scriptFile.name))
condorSubmit.write("Universe            =  vanilla\n")
condorSubmit.write("Requirements        =  Arch==\"X86_64\"  &&  (Machine  !=  \"zebra01.spa.umn.edu\")  &&  (Machine  !=  \"zebra02.spa.umn.edu\")  &&  (Machine  !=  \"zebra03.spa.umn.edu\")  &&  (Machine  !=  \"zebra04.spa.umn.edu\")\n")
condorSubmit.write("+CondorGroup        =  \"cmsfarm\"\n")
condorSubmit.write("getenv              =  True\n")
condorSubmit.write("Request_Memory      =  4 Gb\n")
condorSubmit.write("Log                 =  %s.log\n"%(outDir))

for file in os.listdir(inputDir):

    temp = file.split("_")
    for i in range(len(temp)):

        if temp[i][:5] == "model":
            model = int(temp[i][5:])
        if temp[i][:7] == "version":
            version = int(temp[i][7:])
    
    # Extract output filename from input filename
    outFilename = "Digi_"+str(file)

    condorSubmit.write("Arguments = %s %s %s\n"%(inputDir+"/"+file,outputDir,outFilename))
    condorSubmit.write("Queue\n")

condorSubmit.close()

os.system("chmod u+rwx %s/runDigiJob.sh"%(outDir))
command = "condor_submit " + condorSubmit.name + "\n"
subprocess.call(command.split())
