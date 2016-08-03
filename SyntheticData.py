#	----------------------------------------------
#	----------------------------------------------
#		SYNTHETIC FAULT-SLIP DATA
#		Version: 1.0
#		Author: Prithvi Thakur
#		Last Modified: 03 August 2016, 2:36 p.m.
#	----------------------------------------------
#	----------------------------------------------

#	Generates dip direction, dip angle, trend and plunge of a lineation given a
#	stress tensor. You can change the sigma1, sigma3 and phi values to generate
#	a stress tensor. 
#	Note that plunge is negative for a reverse fault and positive for a normal fault.
#	OUTPUT is an excel file with .xls format. 

import numpy as np
import xlwt

book = xlwt.Workbook()

sheet1 = book.add_sheet("Sheet1")


#--------------------
#	PREALLOCATION
#--------------------

r_matrix = np.zeros((3,3))		# Rotation matrix
phi_tensor = np.zeros((3,3))	# Reduced tensor in the principal reference frame
Data = np.zeros((20,4))			# Output Data Generated

dip_direction = np.zeros((20))	# Dip direction of the fault plane
dip_angle = np.zeros((20))		# Dip angle of the fault plane
slip_azimuth = np.zeros((20))	# Azimuth of lineation
slip_plunge = np.zeros((20))	# Plunge of lineation


#---------------------------------------------------
#	PRINCIPAL STRESS ORIENTATIONS (AZIMUTH/PLUNGE)
#---------------------------------------------------
"""	Azimuth and plunge of sigma1 and sigma3. Sigma2 is perpendicular to both."""
sigma1_az = np.radians(134)
sigma1_pl = np.radians(50)

#sigma2_az = np.radians(90)
#sigma2_pl = np.radians(0)

sigma3_az = np.radians(236)
sigma3_pl = np.radians(10)

phi = 1 - 0.3		# Stress ellipsiod ratio. HERE PHI IS 0.3


#-------------------------------------------
#	DIRECTION COSINES FROM AZIMUTH/PLUNGE
#-------------------------------------------
""" This function returns the direction cosine as a one dimensional array, with its input as
	azimuth and plunge of a line."""
dir_cosine = lambda az, pl: np.array([np.cos(pl)*np.sin(az), np.cos(pl)*np.cos(az), -np.sin(pl)])


#-----------------------------------------
#	AZIMUTH/PLUNGE FROM DIRECTION COSINE
#-----------------------------------------
"""	This function returns the azimuth and plunge as an array, with its input as 
	the direction cosines of a line."""
az_pl = lambda vector: np.array([np.degrees(np.arctan2(vector[0], vector[1])), 
								np.degrees(np.arcsin(-vector[2]/np.linalg.norm(vector)))])


#-------------------
#	ROTATION MATRIX
#-------------------
"""	Rotation matrix is generated using the azimuth and plunge of the principal
	stress axes."""
r_matrix[0,:] = dir_cosine(sigma3_az, sigma3_pl)
a = r_matrix[0,:]
r_matrix[2,:] = dir_cosine(sigma1_az, sigma1_pl)
c = r_matrix[2,:]
r_matrix[1,:] = np.cross(c, a)

#--------------------------------------------
#	CREATING A KNOWN REDUCED STRESS TENSOR
#--------------------------------------------

phi_tensor = np.array([[1, 0, 0],[0, phi, 0],[0, 0, 0]])		# Reduced stress tensor in the principal frame
tensor = np.transpose(r_matrix).dot(phi_tensor.dot(r_matrix))	# Reduced stress tensor in the geographic frames
print tensor


#---------------------------------------------
#	ORIENTATION OF A PLANE- DIP DIRECTION/DIP
#---------------------------------------------
zz = 0
while zz<20:		# Number of fault-slip datum points

	d_deg = np.random.randint(0,360)	# DIP DIRECTION: randomly initialized between 0-360
	d = np.radians(d_deg)				# Dip direction in radians

	flag = 1

	""" The while loop is because for a uniform distribution of strike
		there are more gentle dips than steep dips in the nature. The dip 
		follows a cosine distribution. """
	while flag == 1:
		temp_dip = np.radians(np.random.random_integers(0,90))
		x = np.random.random()
		if x < np.cos(temp_dip):
			p = temp_dip	#	DIP ANGLE
			flag = 2

	normal_vector = dir_cosine(d, p-(np.pi/2))	# Normal vector of a plane
	stress_vector = tensor.dot(normal_vector)	# Stress Vector: from stress tensor and normal vector
	stress_vector_magnitude = np.linalg.norm(stress_vector) # Magnitude of stress vector
	normal_component = stress_vector.dot(normal_vector)		# Normal component of the stress vector
	slip_vector = stress_vector - normal_vector*normal_component	
	shear_component = np.linalg.norm(slip_vector, 2)	# Shear component of the stress vector


	if normal_component/shear_component >=0.85:	# Byerlee's Law
		
		#print slip_vector
		#print normal_vector
		slip_orientation = az_pl(slip_vector)
		dip_direction[zz] = d
		dip_angle[zz] = p
		slip_azimuth[zz] = round(slip_orientation[0] if slip_orientation[0] > 0 else (slip_orientation[0] +360)) 
		slip_plunge[zz] = round(slip_orientation[1])
		zz+=1
		

Data[:,0] = np.degrees(dip_direction)
Data[:,1] = np.degrees(dip_angle)
Data[:,2] = slip_azimuth
Data[:,3] = slip_plunge

for row, array in enumerate(Data):
	for col, value in enumerate(array):
		sheet1.write(row, col, value)
name = "demo.xls"
book.save(name)
print Data
