 #!/usr/bin/env python

import argparse
from time import strftime
import os

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("-f", "--lhefile"    , dest="lhefile"   , help="name of files in directory", required=True)
parser.add_argument("-x", "--xcoordinate", dest="x"         , help="x-coordinate"              , default=0  , required=True, type=float)
parser.add_argument("-y", "--ycoordinate", dest="y"         , help="y-coordinate"              , default=0  , required=True, type=float)
parser.add_argument("-z", "--zcoordinate", dest="z"         , help="z-coordinate"              , default=-50, required=True, type=float)
parser.add_argument("-o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd())
arg = parser.parse_args()

outDir = arg.outDir
nevts = 0
#Check for trailing slash on ouput dir and delete
if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
infileType = "."+arg.lhefile.split(".")[-1]

if not os.path.isdir(outDir):
     print "Output directory does not exist!"
     quit()

if not os.path.isfile(arg.lhefile):
     print "Input file does not exist!"
     quit()

if infileType != ".lhe":
     print "Input file is of incorrect type \"%s\"!"%(infileType)
     quit()

f = open(arg.lhefile, "r")
linesList = f.readlines()
outFilename = arg.lhefile.split("/")[-1].split(".lhe")[0]

lhefile = open("%s/%s_x%s_y%s_z%s_modified.lhe"%(outDir,outFilename,arg.x,arg.y,arg.z), "w")
lhefile.write("<header>\n")
lhefile.write("This file contains signal electrons.\n")
lhefile.write("Do not edit this file manually.\n")
lhefile.write("File created on %s at %s\n"%(strftime("%Y-%m-%d"),strftime("%H:%M:%S")))
lhefile.write("</header>\n")

for j in range(0,len(linesList)):
    parsedLine = linesList[j].split(" ")
    parsedLine = filter(None, parsedLine)

    if parsedLine[0] == "11" and parsedLine[1] == "1":

        px = float(parsedLine[6]) # In GeV
        py = float(parsedLine[7]) # In GeV
        pz = float(parsedLine[8]) # In GeV

        pTot = (px**2 + py**2 + pz**2)**(0.5)

        px_u = px/pTot
        py_u = py/pTot
        pz_u = pz/pTot

        eEnergy = pTot

        KineticE = eEnergy - 0.000510999  

        lhefile.write("<event>\n")
        lhefile.write("11 %g %g %g %g %g %g %g\n"%(arg.x,arg.y,arg.z,px_u,py_u,pz_u,KineticE))
        lhefile.write("</event>\n")

        nevts += 1

lhefile.close()

os.system("mv %s/%s_x%s_y%s_z%s_modified.lhe %s/%s_%s_x%s_y%s_z%s_modified.lhe"%(outDir,outFilename,arg.x,arg.y,arg.z,outDir,nevts,outFilename,arg.x,arg.y,arg.z))

f.close()
