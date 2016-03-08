#			----------------------------------
#			----------------------------------
#				SYNTHETIC FAULT-SLIP DATA
#			----------------------------------
#			----------------------------------

import numpy as np
import xlwt

book = xlwt.Workbook()

sheet1 = book.add_sheet("Sheet1")

#	PREALLOCATION
#--------------------
r_matrix = np.zeros((3,3))
phi_tensor = np.zeros((3,3))
Data = np.zeros((20,4))

dip_direction = np.zeros((20))
dip_angle = np.zeros((20))
slip_azimuth = np.zeros((20))
slip_plunge = np.zeros((20))

#	PRINCIPAL STRESS ORIENTATIONS (AZIMUTH/PLUNGE)
#---------------------------------
sigma1_az = np.radians(134)
sigma1_pl = np.radians(50)

#sigma2_az = np.radians(90)
#sigma2_pl = np.radians(0)

sigma3_az = np.radians(236)
sigma3_pl = np.radians(10)

phi = 1 - 0.3	#	HERE PHI IS 0.3 

#	DIRECTION COSINES FROM AZIMUTH/PLUNGE
#-------------------------------------------
dir_cosine = lambda az, pl: np.array([np.cos(pl)*np.sin(az), np.cos(pl)*np.cos(az), -np.sin(pl)])


#	AZIMUTH/PLUNGE FROM DIRECTION COSINE
#-----------------------------------------
az_pl = lambda vector: np.array([np.degrees(np.arctan2(vector[0], vector[1])), 
								np.degrees(np.arcsin(-vector[2]/np.linalg.norm(vector)))])


#	ROTATION MATRIX
#-------------------
r_matrix[0,:] = dir_cosine(sigma3_az, sigma3_pl)
a = r_matrix[0,:]
r_matrix[2,:] = dir_cosine(sigma1_az, sigma1_pl)
c = r_matrix[2,:]
r_matrix[1,:] = np.cross(c, a)

#	CREATING A KNOWN REDUCED STRESS TENSOR
#--------------------------------------------
phi_tensor = np.array([[1, 0, 0],[0, phi, 0],[0, 0, 0]])
#tensor = r_matrix.dot(phi_tensor.dot(np.transpose(r_matrix)))
tensor = np.transpose(r_matrix).dot(phi_tensor.dot(r_matrix))
print tensor


#	ORIENTATION OF A PLANE- DIP DIRECTION/DIP
#--------------------------------------------
zz = 0
while zz<20:
	d_deg = np.random.randint(0,360) #	DIP DIRECTION
	d = np.radians(d_deg)

	flag = 1

	while flag == 1:
		temp_dip = np.radians(np.random.random_integers(0,90))
		x = np.random.random()
		if x < np.cos(temp_dip):
			p = temp_dip	#	DIP ANGLE
			flag = 2

	normal_vector = dir_cosine(d, p-(np.pi/2))
	stress_vector = tensor.dot(normal_vector)
	stress_vector_magnitude = np.linalg.norm(stress_vector)
	normal_component = stress_vector.dot(normal_vector)
	slip_vector = stress_vector - normal_vector*normal_component
	shear_component = np.linalg.norm(slip_vector, 2)


	if normal_component/shear_component >=0.85:
		
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
name = "uselessdata.xls"
book.save(name)
#book.save(TemporaryFile())
#np.savetxt("SynData1.xls", Data, delimiter="	")
print Data
