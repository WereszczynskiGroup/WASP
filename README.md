# WASP
Writhe Analysis Software Package (WASP)
WASP is a software package for characterizing the topological development of ribbon structures particularly targeted at
at use with molecular dynamics (MD) trajectories.

## Dependencies/Requirements
### General:
Python 2.7 or Python 3+ (recommended)

C++

Bash

### For Amber Trajectories:
Amber/AmberTools, specifically cpptraj

## Configuration
Prior to usage, several actions must be taken to build/configure WASP. 

In order to compile the writhe code, it is necessary to execute the **build.sh**executable in the **writhe\_code** directory which can be done by running the command `sh build.sh` in the **writhe\_code** directory.

The **config.txt** file in the **trajectory\_preprocessing** directory can be modified if the writhe code is planned to be stored somewhere other than its default location or if a large number of trajectories will be analyzed in a different directory. This file is whitespace sensitive and its format should not be altered. 

If the writhe code will be stored in a different directory, simply change the default path for the writhe code to the full path to the new location.

If a large number of trajectories will be analyzed in a different directory, the working directory path can be changed to the directory containing the trajectories in oder to avoid typing the full filepath to each trajectory when running analyses. For example, if `test1.mdcrd, test2.mdcrd` and `test3.mdcrd` are in the directory `/home/temp/trajectories/`, making the working directory `/home/temp/trajectories/` will allow the user to run WASP without having to type `/home/temp/trajectories/` in front of each file name (i.e. not having to type the full file path each time).

In order to analyze Amber mdcrd trajectories, the bash executable **build.sh** in the **trajectory\_preprocessing** directory must be be executed in order to make the **amber_mdcrd_strip.sh** file fully executable by WASP.

## Usage

### User Interface
WASP operates primarily through a CLI user-interface managed by the file **master\_UI.py** in the **trajectory\_preprocessing** directory.

A comprehensive list of arguments can be obtained by running **master\_UI.py** with the arguments **-h** or **--help**, i.e. `python master_UI.py -h`.  

Most arguments available to the user are sufficiently explained to the user by rrunning `python master_UI.py -h`. However, several arguments are worthy of further explanation:

**deletepoints:** This argument specifies how many points should be deleted *on each side* of an axis curve and should be used (carefully) with open/linear DNA structures, i.e. DNA structures whose ends are not joined. For example, specifying `20` for the deletepoints argument will delete 20 points from each end of the axis, deleting a total of 40 points. ***deletepoints should be set to 0 if the DNA structure is closed (circular/knotted) or if the autodelete argument is being used instead***.

**autodelete (-ad):** An alternative to **deletepoints** that attempts to automatically remove stray axis points by starting the axis routine one full helical turn into the DNA structure on both sides. This argument is meant to be used with open/linear DNA structures and **deletepoints** should be set to 0 if **autodelete** is used.

**closed (-c, --closed):** This argument can be used to instruct WASP to use a separate routine for evaluating the polar writhe when analyzing closed DNA structures. At this time, the closed writhe routine and the regular polar writhe routine **(-pw)** are equivalent, so this argument is not necessary. 


