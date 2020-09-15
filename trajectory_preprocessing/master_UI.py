import numpy as np
import math
import sys
import os
import argparse
import subprocess
import shlex
from read_functions import *
from axis import *

#declaring arguments
parser = argparse.ArgumentParser()
parser.add_argument("filetype", choices = ["pdb","mdcrd","oxdna", "general"], help = "trajectory filetype")
parser.add_argument("incrd", help = "input coordinate file name")
parser.add_argument("nbp", help = "number of basepairs in the DNA structure", type = int)
parser.add_argument("startframe", help = "desired trajectory frame to begin analysis on (1 indexed)", type = int)
parser.add_argument("endframe", help = "desired trajectory frame to end analysis on (1 indexed)", type = int)
parser.add_argument("outfile", help = "specifies the names of the output files")
parser.add_argument("--prmtop", help = "specifies name of the topology file (required for amber mdcrd trajectories only)")
parser.add_argument("deleteatoms", help = "specifies the number of atoms to delete from each side of the generated axis curve before analysis (not required for circular/closed DNA structures)", type = int)
parser.add_argument("-s", "--stride", help = "allows for only frames selected by the stride to be analyzed (default = 1, must be greater than or equal to 1)", default = "1", type = int)
parser.add_argument("-d", "--debug", nargs='?', default = False, const = True, help = "writes atoms read from trajectory file and axis curve to files")
parser.add_argument("-x", "--xcol", help = "specifies the column number (zero indexed) corresponding to the x coordinates in the PDB file (required for reading PDB files)", type = int)
parser.add_argument("-pw", "--polarwrithe", nargs='?', default = False, const = True, help = "option to calculate polar writhe")
parser.add_argument("-pws", "--polarwrithestar", nargs='?', default = False, const = True, help = "option to calculate extended polar writhe (Wp*)")
parser.add_argument("-di", "--doubleintegral", nargs='?', default = False, const = True, help = "option to calculate the standard gaussian integral writhe")
parser.add_argument("--smooth", nargs='?', default = False, const = True, help = "use smoothing routine on input axis curve for writhe calculation (recommended), does nothing when used with --doubleintegral")
parser.add_argument("-ad", "--autodelete", nargs='?', default = False, const = True, help = "enable automatic endpoint handling for first frame (do not use for closed structures")
parser.add_argument("-c", "--closed", nargs='?', default = False, const = True, help = "analyze closed curve using -pw or -di")
args = parser.parse_args()

#Checks to see if arguments with """"""optional argument""""" dependencies are satisfied
if args.filetype == "mdcrd" and args.prmtop is None:
	parser.error("The option mdcrd requires an Amber topology file to be specified (--prmtop)")

if args.filetype == "pdb" and args.xcol is None:
	parser.error("The option pdb requires --xcol to be specified")

#main variable initialization
incrd = args.incrd
startframe = args.startframe
endframe = args.endframe
nbp = args.nbp
outfile = args.outfile
deleteatoms = args.deleteatoms
stride = args.stride
xcol = args.xcol
prmtop = args.prmtop
auto_delete = args.autodelete

#Reads file path to the writhe code from the configuration file and strips trailing whitespace
f = open("config.txt")
config_info = f.readlines()
writheCodeFilePath = config_info[1].rstrip()
workingDirectory = config_info[3].rstrip()
#Input formatting check
if writheCodeFilePath[-1] != ("/"):
	writheCodeFilePath += "/"

if workingDirectory[-1] != ("/"):
	workingDirectory += "/"

print("reading " + incrd)

if args.filetype == ("general"):
	subprocess.call(["cp", incrd, writheCodeFilePath + "writhe_input_axis"])
	nframes = (endframe - startframe) + 1

if args.filetype == ("oxdna"):
	x,y,z = read_oxDNA(incrd, startframe, endframe, nbp, stride)
	nframes = len(x)/(2*nbp)
	nframes = int(nframes)

if args.filetype == ("pdb"):
	x,y,z = read_pdb(incrd, startframe, endframe, nbp, stride, xcol)
	nframes = len(x)/(2*nbp)
	nframes = int(nframes)

if args.filetype == ("mdcrd"):

	#runs CPPTRAJ script that strips all non-C1' atoms from the Amber trajectory 
	subprocess.call(shlex.split("%samber_mdcrd_strip.sh %s %s %s %s %s %s"%(str(workingDirectory), str(prmtop), str(incrd), str(nbp), str(startframe), str(endframe), str(stride))))
	
#calculates number of frames in stripped trajectory if the total number of frames is divisible by the stride
#i.e. 1000 frames, stride of 10
#the +1 is to include the boundary frame into the count since cpptraj includes both the start and end frames in the trajectory
#i.e. trajin trajectory.mdrcrd 1 100 1 reads every frame from 1 to 100 including frames 1 and 100, so this would give 100 frames
#note that cpptraj is 1 indexed in terms of frame count

	if ((endframe - startframe+1)%stride == 0):
		nframes = (endframe - startframe+1)/stride

#calculates number of frames in stripped trajectory if the total number of frames is not divisible by the stride (double check this....)
#the  +1 at the end of the line is to include the frame that cpptraj adds to the trajectory if a full stride cannot be completed
#i.e trajin trajectory.mdcrd 1 102 10 reads 11 frames not 10

	else:
		nframes = (endframe - startframe+1)/stride + 1
	
	#explicitly convert nframes to integer
	nframes = int(nframes)

	x,y,z = amber_read("stripped_C1.mdcrd", nbp, nframes)

print("done reading")

if args.filetype != ("general"):
	#convert coordinates to float
	x = x.astype(float)
	y = y.astype(float)
	z = z.astype(float)

	#run axis script
	print("generating axis")
	xyz_coords = read_file(x, y, z, nframes, nbp)
	midpt_ri4, midpt_2 = midpts(xyz_coords, nbp, nframes)
	midpt_axis = midpoint_axis(nbp, nframes, midpt_ri4)
	twist = Twist(nframes, midpt_axis, nbp, xyz_coords)
	wrline_axis, deleteatoms = axis_generate(nbp, nframes, midpt_ri4, twist, deleteatoms, auto_delete)

	#if autodelete is not used, the number of atoms deleted on each side of the helix is = to deleteatoms, so must be
	#multiplied by two to account for both ends of the helix
	if args.autodelete == False:
		deleteatoms = 2*deleteatoms

	f = open("%s"%writheCodeFilePath + "writhe_input_axis","w")

	#writes axis coordinates to file for use with the writhe scripts
	for i in range (nframes):
		for j in range (nbp - deleteatoms):
			f.write("%6f	%6f    %6f\n" % (wrline_axis[i][j][0],wrline_axis[i][j][1],wrline_axis[i][j][2]))


	f.close()
	print("done generating axis")

	#writes axis curve and atoms read from original trajectory file out to xyz files for debugging
	if args.debug == (True):
		print("writing debug")
		f = open (incrd + "debug_backbone.xyz","w")
		for i in range (len(x)/(nbp*2)):
			f.write(str(nbp*2))
			f.write("\n \n")
			for j in range (nbp*2):
				f.write("{0}\t{1}\t{2}\t{3}\n".format('H',x[nbp*i+j],y[nbp*i+j],z[nbp*i+j]))
		f.close()

		f = open(incrd + "_debug_axis.xyz","w")
		for i in range (nframes):
			f.write(str(nbp - deleteatoms))
			f.write("\n \n")
		
			for j in range (nbp - deleteatoms):
	
				f.write("{0}\t{1}\t{2}\t{3}\n".format('H',wrline_axis[i][j][0],wrline_axis[i][j][1],wrline_axis[i][j][2]))
		f.close()
		print("done writing")

#calculates polar writhe (open curve)
if args.polarwrithe == (True) and args.closed == (False):
	print("calculating Wp")
	#checks if smoothing routine is requested
	if args.smooth == (True):
		subprocess.call(shlex.split("%spolarWritheGenTrajectory %swrithe_input_axis %s %s %s %s smooth"%(str(writheCodeFilePath), str(writheCodeFilePath),str(nframes), str(nbp-2*deleteatoms),  str(1), str(outfile) + ".pw")))
	
	else:
		subprocess.call(shlex.split("%spolarWritheGenTrajectory %swrithe_input_axis %s %s %s %s"%(str(writheCodeFilePath), str(writheCodeFilePath),str(nframes), str(nbp-2*deleteatoms), str(1), str(outfile) + ".pw")))


if args.polarwrithestar == (True):
	print("calculating Wp*")
	#checks if smoothing routine is requested
	if args.smooth == (True):
		subprocess.call(shlex.split("%spolarWritheGenTrajectoryStar %swrithe_input_axis %s %s %s %s smooth"%(str(writheCodeFilePath), str(writheCodeFilePath),str(nframes), str(nbp-2*deleteatoms),  str(1), str(outfile) + ".pws")))

	else:
		subprocess.call(shlex.split("%spolarWritheGenTrajectoryStar %swrithe_input_axis %s %s %s %s"%(str(writheCodeFilePath), str(writheCodeFilePath),str(nframes), str(nbp-2*deleteatoms), str(1), str(outfile) + ".pws")))

if args.doubleintegral == (True) and args.closed == (False):
	print("calculating Wr")
	subprocess.call(shlex.split("%sDIWritheGenTrajectory %swrithe_input_axis %s %s %s %s"%(str(writheCodeFilePath), str(writheCodeFilePath),str(nframes), str(nbp-2*deleteatoms), str(1), str(outfile) + ".di")))


#calculates polar writhe (closed curve)
if args.polarwrithe == (True) and args.closed == (True):
	print("calculating Wp")
	#checks if smoothing routine is requested
	if args.smooth == (True):
		subprocess.call(shlex.split("%swritheGenTrajectoryClosed %swrithe_input_axis %s %s %s %s smooth"%(str(writheCodeFilePath), str(writheCodeFilePath),str(nframes), str(nbp-2*deleteatoms),  str(1), str(outfile) + ".pw")))
	
	else:
		subprocess.call(shlex.split("%spolarWritheGenTrajectory %swrithe_input_axis %s %s %s %s"%(str(writheCodeFilePath), str(writheCodeFilePath),str(nframes), str(nbp-2*deleteatoms), str(1), str(outfile) + ".pw")))

#calculates double integral (closed curve)
if args.doubleintegral == (True) and args.closed == (True):
	print("calculating Wr")
	subprocess.call(shlex.split("%sDIWritheGenTrajectoryClosed %swrithe_input_axis %s %s %s %s"%(str(writheCodeFilePath), str(writheCodeFilePath),str(nframes), str(nbp-2*deleteatoms), str(1), str(outfile) + ".di")))
