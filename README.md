# WASP
Writhe Analysis Software Package (WASP)
WASP is a software package for characterizing the topological development of ribbon structures particularly targeted
at (but not limited to) use with molecular dynamics (MD) trajectories.

WASP can be used to calculate the polar writhe (Wp), the extended polar writhe (Wp\*) and the double integral writhe (Wr) of various ribbon structures such as DNA helices.

- [WASP](#wasp)
  * [Dependencies/Requirements](#dependencies-requirements)
    + [General:](#general-)
    + [For Amber Trajectories:](#for-amber-trajectories-)
  * [Configuration](#configuration)
  * [Input File Formats](#input-file-formats)
  * [Usage](#usage)
    + [User Interface](#user-interface)
    + [Examples](#examples)
      - [Analyzing a DNA Minicircle (Amber Trajectory)](#analyzing-a-dna-minicircle--amber-trajectory-)
      - [Analyzing a DNA Minicircle (PDB Trajectory)](#analyzing-a-dna-minicircle--pdb-trajectory-)
      - [Analyzing a DNA Minicircle (General Format Trajectory)](#analyzing-a-dna-minicircle--general-format-trajectory-)
      - [Analyzing a linear DNA Plectoneme (oxDNA Trajectory, autodelete)](#analyzing-a-linear-dna-plectoneme--oxdna-trajectory--autodelete-)
      - [Analyzing a linear DNA Plectoneme (oxDNA Trajectory, deleteatoms)](#analyzing-a-linear-dna-plectoneme--oxdna-trajectory--deleteatoms-)
  * [Output Formats](#output-formats)
  * [References](# References)

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

## Input File Formats

Currently, WASP is equipped to handle four different file formats:

1. Amber mdcrd trajectories (with associated topology (prmtop) files)
2. oxDNA xyz trajectories ***Note: WASP will look for any and all atoms named C in the trajectory file to calculate the axis curve. If no C atoms are present, it will not calculate an axis curve. Similarly, if any other structures are present that contain C atoms that you are not trying to analyze, it may try to use those also***
3. PDB trajectories ***Note: WASP will look for any and all atoms named C1' in the trajectory file to calculate the axis curve. If no C atoms are present, it will not calculate an axis curve. Similarly, if any other structures are present that contain C atoms that you are not trying to analyze, it may try to use those also***
4. "General format" trajectories (see below)

If you have a molecular dynamics trajectory that is not one of these formats, the trajectory can be made "WASP friendly" by converting the trajectory to PDB format using VMD. 

Examples of each format can be found in the **examples** directory for reference.

Using the general format option tells WASP that the trajectory you are giving it is ***already an axis curve***. In other words, you are feeding WASP the exact curve that you want to calculate the writhe of, and the axis curve calculation is skipped. This allows WASP to be used more generally if you wish to calculate an axis curve for your DNA structure using some other method prior to using WASP to calculate the writhe or if you want to analyze the writhing of a curve that is not DNA or not a result of a MD trajectory (i.e. coronal flux ropes or rod filaments). This format consists of a tab spaced cartersian (x,y,z) coordinates on each line with no whitespace between lines. It is important to note that the number of points on the curve should be the same in every frame, as WASP has no way of telling otherwise. When using general format trajectories, some functionalities of WASP are unable to be used such as the ability to stride. See the general format example below for more information.


For PDB files, the --xcol argument must be specified. The user should provide the column number (zero-indexed) in the PDB file that contains the x coordinates. In the example PDB trajectory in the **examples** directory, the --xcol argument should be 6.

## Usage

### User Interface
WASP operates primarily through a CLI user-interface managed by the file **master\_UI.py** in the **trajectory\_preprocessing** directory.

A comprehensive list of arguments can be obtained by running **master\_UI.py** with the arguments **-h** or **--help**, i.e. `python master_UI.py -h` or by running `python master_UI.py` with no arguments at all.

Most arguments available to the user are sufficiently explained to the user by rrunning `python master_UI.py -h`. However, several arguments are worthy of further explanation:

**deleteatoms:** This argument specifies how many points should be deleted *on each side* of an axis curve and should be used (carefully) with open/linear DNA structures, i.e. DNA structures whose ends are not joined. For example, specifying `20` for the deleteatoms argument will delete 20 points from each end of the axis, deleting a total of 40 points. ***deleteatoms should be set to 0 if the DNA structure is closed (circular/knotted) or if the autodelete argument is being used instead***.

**autodelete (-ad):** An alternative to **deleteatoms** that attempts to automatically remove stray axis points by starting the axis routine one full helical turn into the DNA structure on both sides. This argument is meant to be used with open/linear DNA structures and **deleteatoms** should be set to 0 if **autodelete** is used.

**closed (-c, --closed):** This argument can be used to instruct WASP to use a separate routine for evaluating the polar writhe when analyzing closed DNA structures. At this time, the closed writhe routine and the regular polar writhe routine **(-pw)** are equivalent, so this argument is not necessary. 

**smooth (--smooth):** Enables a smoothing routine used to improve accuracy when calculating Wp and Wp\* (recommended)

**debug (-d, --debug):** Writes the backbone atoms used for axis curve calculation and the generated axis curve to separate trajectory files. These trajectory files are XYZ trajectory files that meet the specifications required to be read by VMD. This option is useful for determining how many atoms to delete from an axis curve if analyzing a linear DNA helix with the **deleteatoms** option.

### Examples

Sample trajectory files that can be used to test WASP can be found in the **examples** directory. We will demonstrate how to use the functionalities of WASP using these example trajectories. For these examples, we assume that the config.txt file is left as the default.

#### Analyzing a DNA Minicircle (Amber Trajectory)

For the first example, we will analyze a DNA minicircle trajectory genereated from an Amber MD simulation. The file **minicircle.mdcrd** in the **examples** directory represents a 100 frame trajectory of a 108 base pair minicircle. We will accomplish the analysis by running (in the trajectory\_preprocessing directory):

`python master_UI.py mdcrd ../examples/minicircle.mdcrd 108 30 100 minicircle_data 0 --prmtop ../examples/minicircle.prmtop -pw -di --smooth -c -d`

Running this will result in the following:

1. The **minicircle.mdcrd** trajectory will be read into cpptraj along with the **minicircle.prmtop** topology file
2. All atoms will be stripped from the trajectory except for the C1' atoms on the DNA backbone necessary for calculating the DNA axis
3. Frames [30, 100] will be analyzed for a total of 71 frames starting on frame 30 and ending on frame 100
4. The C1' atoms on the backbone and the axis curve will both be written to trajectory files named **minicircle.mdcrd\_debug\_backbone.xyz** and **minicircle.mdcrd\_debug\_axis.xyz** respectively in the **examples** directory
5. No atoms will be deleted from the axis curve prior to writhe analysis
6. The polar writhe (Wp) and the double integral writhe will be calculated for each of the 70 frames and written to files named **minicircle\_data.pw** and **minicircle\_data.di** respectively
7. The smoothing routine will be used for the Wp calculation
8. The closed Wp routine will be used for the Wp calculation

#### Analyzing a DNA Minicircle (PDB Trajectory)

For this example, we will analyze a DNA minicircle trajectory stored in PDB format. This trajectory is identical to the trajectory in the previous (Amber Trajectory) example with 100 frams and 108 base pairs except for the fact that it is in PDB format. This time, we run (in the trajectory\_preprocessing directory):

`python master_UI.py pdb ../examples/minicircle.pdb 108 1 100 minicircle_data 0 -pw -di --smooth -c --xcol 6 --stride 10 -d`

Running this will result in the following:

1. The C1' atoms on the minicircle backbone will be read in every 10th frame (**--stride 10**) between frames [1,100] and used to calculate an axis curve
2. WASP will look for the x-coordinates in the 6th column (zero indexed) of the PDB file
2. Every 10th frame on the interval [1, 100] will be analyzed for a total of 10 frames
3. The C1' atoms on the backbone and the axis curve will both be written to trajectory files named **minicircle.pdb\_debug\_backbone.xyz** and **minicircle.pdb\_debug\_axis.xyz** respectively in the **examples** directory
4. No atoms will be deleted from the axis curve prior to writhe analysis
5. The polar writhe (Wp) and the double integral writhe will be calculated for each of the 10 frames and written to files named **minicircle\_data.pw** and **minicircle\_data.di** respectively
6. The smoothing routine will be used for the Wp calculation
7. The closed Wp routine will be used for the Wp calculation

#### Analyzing a DNA Minicircle (General Format Trajectory)

For this example, we will analyze a DNA minicircle ***axis curve trajectory*** generated from the same minicircle trajectory as used in the previous examples. Since the general format assumes that the file give is already an axis curve, WASP will skip straight to the writhe caluculation routines and many of the WASP functionalities are not able to be used. For general trajectories, the **nbp** argument is simply the number of points on the axis curve in each frame. This time, we run (in the trajectory\_preprocessing directory):

`python master_UI.py general ../examples/minicircle.general 108 1 100 minicircle_data 0 -pw -di --smooth -c`

Running this will result in the following:

1. The entire trajectory will be read in from frame 1 to 100 ([1, 100]). **startframe** and **endframe** must ***always*** be the very first and very last frames in a general format trajectory.
2. No atoms will be deleted from the axis curve prior to writhe analysis
3. The polar writhe (Wp) and the double integral writhe will be calculated for each frame and written to files named **minicircle\_data.pw** and **minicircle\_data.di** respectively
4. The smoothing routine will be used for the Wp calculation
5. The closed Wp routine will be used for the Wp calculation

#### Analyzing a linear DNA Plectoneme (oxDNA Trajectory, autodelete)

For this example, we will analyze a linear DNA plectoneme trajectory generated from an oxDNA simulation. This trajectory has a total of 120 frames and the DNA helix has 600 basepairs. Because it is not a closed DNA structure like a minicircle, we must use a routine (**autodelete** in this example) to remove stray axis points generated by the WrLINE algorithm. This time, we run (in the trajectory\_preprocessing directory):

`python master_UI.py oxdna ../examples/plectoneme.oxdna 600 1 120 plectoneme_data 0 -pw -pws --smooth -ad -d --stride 10`

Running this will result in the following:

1. The C (not C1', this is an oxDNA convention) atoms on the backbone will be read every 10th frame (**--stride 10**) between frames [1,120] and used to calculate an axis curve
2. The axis curve routine will calculate an axis curve for each frame in the stride for a total of 12 frames
3. The autodelete routine (**-ad**) will attempt to remove all stray axis points from the calculated axis curves
4. Every 10th frame on the interval [1, 12] will be analyzed for a total of 12 frames
5. The C atoms on the backbone and the axis curve will both be written to trajectory files named **plectoneme.oxdna\_debug\_backbone.xyz** and **plectoneme.oxdna\_debug\_axis.xyz** respectively in the **examples** directory
6. The polar writhe (Wp) and the extended polar writhe (Wp\*) will be calculated fo each frame included in the stride and written to files named **plectoneme\_data.pw** and **plectoneme_data.pws** respectively
7. The smoothing routine will be used for the Wp and Wp\* calculations

#### Analyzing a linear DNA Plectoneme (oxDNA Trajectory, deleteatoms)

For this example, we will analyze a linear DNA plectoneme trajectory generated from an oxDNA simulation. This trajectory has a total of 120 frames and the DNA helix has 600 basepairs. Because it is not a closed DNA structure like a minicircle, we must use a routine (**deleteatoms** in this example) to remove stray axis points generated by the WrLINE algorithm. This time, we run (in the trajectory\_preprocessing directory):

`python master_UI.py oxdna ../examples/plectoneme.oxdna 600 1 120 plectoneme_data 25 -pw -pws --smooth -d --stride 10`

Running this will result in the following:

1. The C (not C1', this is an oxDNA convention) atoms on the backbone will be read every 10th frame (**--stride 10**) between frames [1,120] and used to calculate an axis curve
2. The axis curve routine will calculate an axis curve for each frame in the stride for a total of 12 frames
3. 25 atoms will be deleted from each end of the axis curve for a total of 50 atoms deleted to remove stray axis points (**deleteatoms = 25**)
4. Every 10th frame on the interval [1, 12] will be analyzed for a total of 12 frames
5. The C atoms on the backbone and the axis curve will both be written to trajectory files named **plectoneme.oxdna\_debug\_backbone.xyz** and **plectoneme.oxdna\_debug\_axis.xyz** respectively in the **examples** directory
6. The polar writhe (Wp) and the extended polar writhe (Wp\*) will be calculated fo each frame included in the stride and written to files named **plectoneme\_data.pw** and **plectoneme_data.pws** respectively
7. The smoothing routine will be used for the Wp and Wp\* calculations

## Output Formats

Output files that result from Wp or Wp\* have five components and look like the following:

```
-0.0230357 0.0235573 0.000521591 {0} 1
0.274211 0.160476 0.434687 {0} 2
0.273186 0.247644 0.520831 {0 0 0 0 0 0} 3
0.313301 0.458862 0.772162 {0 0 0 0 0 0} 4
0.319946 0.814396 1.13434 {0 0 0 0 0 0} 5
0.141583 1.06685 1.20843 {0} 6
0.260391 1.13204 1.39243 {0} 7
0.147154 1.22889 1.37604 {0} 8

```

For both Wp and Wp\* output files, the first column corresponds to the local component of the writhe (Wpl), the second column corresponds to the nonlocal component of the writhe (Wpnl) and the third column corresponds to the total Wp or Wp\* value (Wp + Wpnl). The last column indicates the frame number. The numbers between the curly brackets list the integer windings for every subsection pair of the curve, but are likely irrelevant for most analyses.

Output files that result from double integral writhe calculations have 2 components and look like the following:

```
0.00143835 1
0.534908 2
0.604642 3
0.791726 4
1.11398 5
1.28044 6
1.42515 7
1.42484 8

```

In this case, the first column is the writhe value and the second column is the frame number.

## References

WASP uses an implementation of the WrLINE method developed by Thana Sutthibutpong, Sarah A. Harris, and Agnes Noy (https://doi.org/10.1021/acs.jctc.5b00035) for generating axis curves of DNA helices.
