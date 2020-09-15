# WASP
Writhe Analysis Software Package (WASP)
WASP is a software package for characterizing the topological development of ribbon structures particularly targeted at
at use with molecular dynamics (MD) trajectories.

## Dependencies/Requirements
### General:
Python 2.7 or Python 3+ (recommended)
C++
Bash

###For Amber Trajectories:
Amber/AmberTools, specifically cpptraj

## Configuration
Prior to usage, several actions must be taken to build/configure WASP. 

In order to compile the writhe code, it is necessary to execute the **build.sh**executable in the **writhe\_code** directory which can be done by running the command `sh build.sh` in the **writhe\_code** directory.

The **config.txt** file in the **trajectory\_preprocessing** directory can be modified if the writhe code is planned to be stored somewhere other than its default location or if a large number of trajectories will be analyzed in a different directory. 

If the writhe code will be stored in a different directory, simply change the default path for the writhe code to the full path to the new location.

If a large number of trajectories will be analyzed in a different directory, the working directory path can be changed to the directory containing the trajectories in oder to avoid typing the full filepath to each trajectory when running analyses. For example, if test1.mdcrd, test2.mdcrd and test3.mdcrd are in the directory /home/temp/trajectories/, making the working directory /home/temp/trajectories/ will allow the user to run WASP without having to type /home/temp/trajectories/ in front of each file name (i.e. not having to type the full file path each time).

In order to analyze Amber mdcrd trajectories, the bash executable **build.sh** in the **trajectory\_preprocessing** directory must be be executed in order to make the **amber_mdcrd_strip.sh** file fully executable by WASP.

## Usage

WASP operates primarily through a CLI user-interface managed by the file **master\_UI.py** in the **trajectory\_preprocessing** directory.
