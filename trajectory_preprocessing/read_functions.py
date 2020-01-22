import numpy as np

def read_oxDNA(incrd, startframe, endframe, nbp, stride):

	xcrd = []
	ycrd = []
	zcrd = []

	nbp2 = 2*nbp

	#read only C atoms from oxDNA XYZ trajectory
	#WARNING: MAKE SURE THERE ARE NOT MULTIPLE MOLECULES CONTAINING C ATOMS! ONLY A SINGLE DNA STRUCTURE IS SUPPORTED

	f = open(incrd,"r")

	for line in f:
		if line.startswith('C'):
			data_line = line.rstrip().split(' ')
			while ('') in data_line: data_line.remove('')
			xcrd.append(data_line[1])
			ycrd.append(data_line[2])
			zcrd.append(data_line[3])

	xcrd = np.asarray(xcrd)
	ycrd = np.asarray(ycrd)
	zcrd = np.asarray(zcrd)

	#adjusting for zero indexing
	frame_start = startframe - 1
	frame_end = endframe - 1

	#create list of coordinates only including coordinates within the specified frame range
	xframes = xcrd[nbp2*frame_start:nbp2*frame_end+nbp*2]
	yframes = ycrd[nbp2*frame_start:nbp2*frame_end+nbp*2]
	zframes = zcrd[nbp2*frame_start:nbp2*frame_end+nbp*2]

	#handles strides
	if stride != 1:
		xframes_stride = np.array([])
		yframes_stride = np.array([])
		zframes_stride = np.array([])
		for i in range (frame_end - frame_start):

			#add coordinates to output file if they are in a frame that is included in the stride i.e. if stride = 2 and i = 4, include coordinates
		       
			#if stride = 2 and i = 5, don't include coordinates
			if i%stride == 0:
				xframes_stride = np.append(xframes_stride,xframes[i*nbp2:i*nbp2+nbp2])
				yframes_stride = np.append(yframes_stride,yframes[i*nbp2:i*nbp2+nbp2])
				zframes_stride = np.append(zframes_stride,zframes[i*nbp2:i*nbp2+nbp2])

		return xframes_stride,yframes_stride,zframes_stride
	else:

		return xframes, yframes, zframes

def read_pdb(incrd, startframe, endframe, nbp, stride, xcol):

        xcrd = []
        ycrd = []
        zcrd = []

        nbp2 = 2*nbp

        #read only C1' atoms from PDB
        #WARNING: MAKE SURE THERE ARE NOT MULTIPLE MOLECULES CONTAINING C1' ATOMS! ONLY A SINGLE DNA STRUCTURE IS SUPPORTED

        f = open(incrd,"r")

        for line in f:
                if line.startswith('ATOM') and "C1'" in line:
                        data_line = line.rstrip().split(' ')
                        while ('') in data_line: data_line.remove('')
                        xcrd.append(data_line[xcol])
                        ycrd.append(data_line[xcol+1])
                        zcrd.append(data_line[xcol+2])

        xcrd = np.asarray(xcrd)
        ycrd = np.asarray(ycrd)
        zcrd = np.asarray(zcrd)

        #adjusting for zero indexing
        frame_start = startframe - 1
        frame_end = endframe - 1

        #create list of coordinates only including coordinates within the specified frame range
        xframes = xcrd[nbp2*frame_start:nbp2*frame_end+nbp*2]
        yframes = ycrd[nbp2*frame_start:nbp2*frame_end+nbp*2]
        zframes = zcrd[nbp2*frame_start:nbp2*frame_end+nbp*2]

        #handles strides
        if stride != 1:
                xframes_stride = np.array([])
                yframes_stride = np.array([])
                zframes_stride = np.array([])
                for i in range (frame_end - frame_start):

                        #add coordinates to output file if they are in a frame that is included in the stride i.e. if stride = 2 and i = 4, include coordinates
                        #if stride = 2 and i = 5, don't include coordinates
                        if i%stride == 0:
                                xframes_stride = np.append(xframes_stride,xframes[i*nbp2:i*nbp2+nbp2])
                                yframes_stride = np.append(yframes_stride,yframes[i*nbp2:i*nbp2+nbp2])
                                zframes_stride = np.append(zframes_stride,zframes[i*nbp2:i*nbp2+nbp2])

                return xframes_stride,yframes_stride,zframes_stride
        else:

                return xframes, yframes, zframes



def amber_read(incrd, nbp, nframes):

        f = open(incrd, "r")

        #ignore first line containing "Cpptraj Generated Trajectory" header
        f.readline()

        #list of all coordinates before distinguishing between x,y,z    
        allcoords=[]

        #extracting coordinates from mdcrd format
        for line in f:

                line_length = len(line)
                n = int(line_length/8)

                for i in range (0,n):

                        allcoords.append(float(line[i*8:(i+1)*8]))

        #make coordinate list an array for slicing and reshape to get x,y,z coordinates in columns 0,1,2
        allcoords = np.asarray(allcoords)
        allcoords = np.reshape(allcoords,(nframes*nbp*2,3))

        #slice coordinate array to get x,y,z coordinates from columns
        xcrd = allcoords[:,0]
        ycrd = allcoords[:,1]
        zcrd = allcoords[:,2]

        f.close()

        return xcrd, ycrd, zcrd




