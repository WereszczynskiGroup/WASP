import numpy as np
import math
import sys
import os

#Implementation of the WrLINE axis method (Thana Sutthibutpong, Sarah Harris and Agnes Noy  https://doi.org/10.1021/acs.jctc.5b00035)
def read_file(x, y, z, nframes, nbp):
	
	nbp2 = 2*nbp

	xyz_coords = []
	for i in range(0,nframes):
		xyz_frame = []
		for j in range(0,nbp2):
			xyz_frame.append([x[j+nbp2*i],y[j+nbp2*i],z[j+nbp2*i]])
	
		xyz_coords.append(xyz_frame)

	xyz_coords = np.array(xyz_coords)
	return xyz_coords


def midpts(xyz_coords, nbp, nframes):

	nbp2 = 2*nbp-1
	midpt_ri4 = []
	midpt_2 = []
	for i in range (0,nframes):
		midpts_4 = []
		midpts_2 = []
		for j in range (0,nbp):
			midpts_2.append((xyz_coords[i][nbp2-j]+xyz_coords[i][j])/2)	
			
			if j == (nbp-1):
				ri = (xyz_coords[i][nbp2-j]+xyz_coords[i][j]+xyz_coords[i][nbp2]+xyz_coords[i][0])/4			

			else:
				ri = (xyz_coords[i][nbp2-j]+xyz_coords[i][j]+xyz_coords[i][j+1]+xyz_coords[i][nbp2-j-1])/4
		
				
			midpts_4.append(ri)

		midpt_ri4.append(midpts_4)
		midpt_2.append(midpts_2)

	return midpt_ri4, midpt_2



def rotations(r_ab, r_cd):

	r_ab = r_ab/np.linalg.norm(r_ab)

def midpoint_axis(nbp, nframes, midpt_ri4):

	midpt_axis = []

	for i in range (0, nframes):
		
		midpt_axis_frame = []
		#pt_sum = np.zeros(3)
		for j in range (0,nbp):
			pt_sum = np.zeros(3)	
			pt_sum += midpt_ri4[i][j]
		
			for k in range (1,6):
				if (j+k) > (nbp-1) and (j-k) < 0:
					pt_sum += midpt_ri4[i][j+k-nbp]+midpt_ri4[i][j-k+nbp]

				elif (j+k) > (nbp-1) and (j-k) >= 0:
					pt_sum += midpt_ri4[i][j+k-nbp]+midpt_ri4[i][j-k]

				elif (j-k) < 0 and (j+k) <= (nbp-1):

					pt_sum += midpt_ri4[i][j+k]+midpt_ri4[i][j-k+nbp]
				else:
					pt_sum += midpt_ri4[i][j+k]+midpt_ri4[i][j-k]

			#11 = 2x5 + 1 for neigboring basepairs and the inclusive basepair
			ax_pt = pt_sum/(11)
			midpt_axis_frame.append(ax_pt)
		
		midpt_axis.append(midpt_axis_frame)
	
	return midpt_axis			
					
def Twist(nframes, midpt_axis, nbp, xyz_coords):

	nbp2 = 2*nbp-1
	twist = []
	for i in range (0,nframes):
		twist_frame = []
		for j in range (0,nbp):
			if (j+1) > (nbp-1):
				zi = midpt_axis[i][j+1-nbp]-midpt_axis[i][j-1]
				r_ab = xyz_coords[i][nbp2-j] - xyz_coords[i][j]
				r_cd = xyz_coords[i][nbp2] - xyz_coords[i][0]
			else:
				zi = midpt_axis[i][j+1]-midpt_axis[i][j-1]
				r_ab = xyz_coords[i][nbp2-j] - xyz_coords[i][j]
				r_cd = xyz_coords[i][nbp2-j-1] - xyz_coords[i][j+1]

			r_ab = r_ab/np.linalg.norm(r_ab)
			r_cd = r_cd/np.linalg.norm(r_cd)
			
			zi = zi/np.linalg.norm(zi)

			alpha = math.atan2(zi[1],zi[2])
			beta = np.arctan(-1*zi[0]/np.sqrt(zi[1]**2+zi[2]**2))
			rotate_x = np.array([[1,0,0],[0,np.cos(alpha),-1*np.sin(alpha)],[0,np.sin(alpha),np.cos(alpha)]])
			rotate_y = np.array([[np.cos(beta),0,np.sin(beta)], [0,1,0], [-1*np.sin(beta), 0, np.cos(beta)]])

			rotxy = np.dot(rotate_y,rotate_x)
			
			r_ab = np.dot(rotxy, r_ab)
			r_cd = np.dot(rotxy, r_cd)

			gamma = -1*math.atan2(r_ab[1],r_ab[0])
			rotz = np.array([[np.cos(gamma),-1*np.sin(gamma),0], [np.sin(gamma), np.cos(gamma), 0], [0,0,1]])
			
			r_cd = np.dot(rotz, r_cd)
			
			tw_angle = math.atan2(r_cd[1],r_cd[0])*180/np.pi
			twist_frame.append(tw_angle)
		twist.append(twist_frame)
	return twist

def axis_generate(nbp, nframes, midpt_ri4, twist, deleteatoms):
	wrline_axis = []

	for i in range (0,nframes):
		wrline_frame = []
		for j in range (0,nbp):
			theta_m = 0
			pt_sum = np.zeros(3)
			
			theta_m += twist[i][j]
			pt_sum += midpt_ri4[i][j]
			
			k = 0

			while (theta_m < 360.0):
				k += 1
				theta_mminus = theta_m
				
				if (j+k) > (nbp-1):
					if (j-k) < 0:
						print(j-k+nbp,j+k-nbp)
						pt_sum += midpt_ri4[i][j-k+nbp] + midpt_ri4[i][j+k-nbp]
						theta_m += twist[i][j-k+nbp] + twist[i][j+k-nbp]
					else:
						pt_sum += midpt_ri4[i][j-k] + midpt_ri4[i][j+k-nbp]
						theta_m += twist[i][j-k] + twist[i][j+k-nbp]	
				else:
					if (j-k) < 0:
						pt_sum += midpt_ri4[i][j-k+nbp] + midpt_ri4[i][j+k]
						theta_m += twist[i][j-k+nbp] + twist[i][j+k]
					else:
						pt_sum += midpt_ri4[i][j-k] + midpt_ri4[i][j+k]
						theta_m += twist[i][j-k] + twist[i][j+k]
			weight = (360 - theta_mminus)/(theta_m - theta_mminus)
			if (j-k) < 0 and (j+k) > (nbp-1):
				pt_sum = pt_sum - (1-weight)*(midpt_ri4[i][j-k+nbp] + midpt_ri4[i][j+k-nbp])
			elif (j-k) < 0 and (j+k) <= (nbp-1):
				pt_sum = pt_sum - (1-weight)*(midpt_ri4[i][j-k+nbp] + midpt_ri4[i][j+k])

			elif (j+k) > (nbp-1) and (j-k) >= 0:
				pt_sum = pt_sum - (1-weight)*(midpt_ri4[i][j-k] + midpt_ri4[i][j+k-nbp])
			else:
				pt_sum = pt_sum - (1-weight)*(midpt_ri4[i][j-k] + midpt_ri4[i][j+k])

			h_i = pt_sum/(2*(k+weight)-1)
			wrline_frame.append(h_i)
		wrline_axis.append(wrline_frame)

	if deleteatoms != 0:
		for i in range(nframes):
			del wrline_axis[i][nbp-deleteatoms:nbp]
			del wrline_axis[i][0:deleteatoms]

	return wrline_axis


