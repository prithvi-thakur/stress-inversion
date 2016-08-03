#	----------------------------------------------------------
#	----------------------------------------------------------
# 	 	THE GA (GENETIC ALGORITHM) METHOD FOR STRESS INVERSION 
#		FROM HOMOGENEOUS/SINGLE PHASE FAULT-SLIP DATA
#		Version: 1.0
#		Author: Prithvi Thakur
#		Last Modified: 03 August 2016, 2:51 p.m. 
#	----------------------------------------------------------
#	----------------------------------------------------------

#	Using the genetic algorithm, this program inverts the fault-slip data
#	to obtain a reduced stress tensor. Line 71 takes the input data in excel 
#	file format with four columns: dip direction, dip angle, trend and plunge
#	of lineations. The input data file should be in the same directory as this
#	program. You can play with the variables in lines 179-183 to tweak the algorithm
#	and possibly obtain better inversion results.

import numpy as np 
import random
import math
import xlrd
import matplotlib.pyplot as plt

random.seed()	#seed the random number - system time


#-------------------
#	ROTATION MATRIX	 (3-dimensional rotation matrix: check out wikipedia)
#-------------------
""" x, y, z are the euler angles, euler_rotation multiplies the rotations about the three axes 
	to get a conposite rrotation matrix """
euler_z = lambda z: np.array([[np.cos(z), -np.sin(z), 0], [np.sin(z), np.cos(z), 0] ,[0, 0, 1]])
euler_y = lambda y: np.array([[np.cos(y), 0, np.sin(y)], [0, 1, 0], [-np.sin(y), 0, np.cos(y)]])
euler_x = lambda x: np.array([[1, 0, 0], [0, np.cos(x), -np.sin(x)], [0, np.sin(x), np.cos(x)]])

euler_rotation = lambda x, y, z: z.dot(y.dot(x))


#-----------------------------------------
#	DIRECTION COSINES FROM AZIMUTH/PLUNGE
#-----------------------------------------
""" This function returns the direction cosine as a one dimensional array, with its input as
	azimuth and plunge of a line."""
dir_cosine = lambda az, pl: np.array([np.cos(pl)*np.sin(az), np.cos(pl)*np.cos(az), -np.sin(pl)])


#----------------------------------------
#	AZIMUTH/PLUNGE FROM DIRECTION COSINE
#----------------------------------------
"""	This function returns the azimuth and plunge as an array, with its input as 
	the direction cosines of a line."""
az_pl = lambda vector: np.array([np.degrees(np.arctan2(vector[0], vector[1])), 
								np.degrees(np.arcsin(-vector[2]/np.linalg.norm(vector)))])

#--------------------------------------------------
#	DECIMAL TO BINARY FUNCTION, b = number of bits
#--------------------------------------------------
dec2bin = lambda x, n:format(x, 'b').zfill(n) 


#---------------------
#	FITNESS FUNCTION
#---------------------
def fitnessFunc(stress_tensor, pop_size, n_i):
	"""	This function computes the mean of the misfits of the individual data points.  
		Misfit is calculated as the difference between observed and the theoretical
		shear stress directions. For the mathematics of this function refer to my 
		'soon to be published paper' """

	input_data = xlrd.open_workbook('demo.xls')	# change the filename here for the data to be inverted
	first_sheet = input_data.sheet_by_index(0)
	dip_d = np.radians(first_sheet.col_values(0))	# Dip direction
	dip_a = np.radians(first_sheet.col_values(1))	# Dip angle
	trend = np.radians(first_sheet.col_values(2))	# Trend of lineation
	plunge = np.radians(first_sheet.col_values(3))	# Plunge of lineation
	length = len(dip_d)		# Number of data points

	#	PREALLOCATION
	normal_vector = np.zeros((length,3))	# Normal vector of a plane
	slip_vector = np.zeros((length,3))		# Slip vector of the lineation
	fitness = np.zeros((pop_size,1))		# Measure of misfit
	
	#	SUB FITNESS FUNCTION
	def subfitnessFunc(sigma, slip, normal):
		"This function computes the composite misfit for each datum"
		stress_vector = sigma.dot(normal)	# Stress Vector: from stress tensor and normal vector
		observed_shear_component = stress_vector.dot(slip)	# Observed Shear component of stress vector
		stress_vector_mag = np.linalg.norm(stress_vector, 2)	# Stress vector magnitude
		normal_component = stress_vector.dot(normal)		# Normal component of stress vector
		shear_component = stress_vector - (normal * normal_component)	# Theoretically obtained shear component of stress vector
		unit_shear_component = shear_component / np.linalg.norm(shear_component, 2)

		new_misfit = np.absolute(np.arccos(unit_shear_component.dot(slip)))	# Angular misfit between theretical and observed component

		f2 = np.sin((new_misfit/2))	#f2: sine of angular misfit
		f1 = np.absolute( np.absolute(observed_shear_component) - np.sqrt(np.square(stress_vector_mag) - np.square(normal_component)))

		lamb = 0.7	# Weighting parameter
		f = lamb*f1 + (1-lamb)*f2	# Composite Misfit function for one datum point

		return f

	#	COMPUTING NORMAL AND SHEAR VECTOR
	for i in range(0, length):	
		normal_vector[i,:] = dir_cosine(dip_d[i], dip_a[i]-(np.pi/2))
		slip_vector[i,:] = dir_cosine(trend[i], plunge[i])

	#	COMPUTING FITNESS ARRAY
	for q in range(0,pop_size):
		sum_f = 0

		for w in range(0,length):
			normal = normal_vector[w,:]
			slip = slip_vector[w,:]
			sigma = stress_tensor[q,:,:]	

			f = subfitnessFunc(sigma, slip, normal)
			sum_f += f 
		
		fitness[q] = sum_f/length	# mean of composite misfit functionS

	return fitness
#------------------	


#---------------------
#	CROSSVER FUNCTION
#---------------------
def crossoverFunc(parents, size, bits):
	"""This function performs single point crossover for input array of strings.
	Sixty percent individuals undergo crossover"""
	children = np.zeros((size), np.dtype('a6'))

	for i in range(0, int(size/2)):
		x_site = np.random.randint(0, bits - 1)	# Crossover site: between 0 and the number of bits
		x1 = parents[i]		# First parent for crossover
		x2 = parents[size - i - 1]	# Second parent for crossover

		if (np.random.randint(0, 100)) > 40 :	# Crossover probability
			ch1 = x1[0:x_site] + x2[x_site:bits]	# Exchanging the bits
			ch2 = x2[0:x_site] + x1[x_site:bits]

			children[i] = ch1
			children[size - i - 1] = ch2
		
		else:		# Children not undergoing crossover
			children[i] = x1
			children[size - i - 1] = x2


	return children
#-------------------


#---------------------
#	MUTATION FUNCTION
#---------------------
def mutationFunc(parents, size, bits):
	"This function performs mutation for input array of strings"
	child = np.dtype('a6')
	prob = np.random.randint(0,100)

	if prob > 95:		# Probability of mutation
		s = np.random.randint(0, size)	# Parent undergoing mutation
		b = np.random.randint(0, bits)	# Bit of that parent undergoing change
		p = parents[s]
		temp1 = p[0:b] 
		temp2 = str(1 - int(p[b]))
		temp3 = p[b+1:bits]
		child = temp1 + temp2 + temp3
		parents[s] = child

	return parents


#-----------------
#   PREALLOCATION
#-----------------
pop_size = 1000		# Population size
noi = 60			# Number of Iterations
bits = 6			# Number of bits of each parameter
elite = 5			# Percentage of elite individuals
t = int(elite * 0.01 * pop_size)	# Fraction of elite individuals 
print "The Population size is %d and the number of iterations is %d" % (pop_size, noi)

reduced_tensor = np.zeros((pop_size,3,3))	# Reduced tensor in principal frame: with only phi in it
stress_tensor = np.zeros((pop_size,3,3))	# Reduced tensor in geographical frame: phi and the euler angles
t_matrix = np.zeros((pop_size,3,3))			# Transformation matrix / Rotation matrix

fitness = np.zeros((pop_size))		# Array of misfits

m_alpha = np.zeros((pop_size - t))	#----- 
m_beta = np.zeros((pop_size - t))	# Mating pool
m_gamma = np.zeros((pop_size - t))	#
m_phi = np.zeros((pop_size - t))	#----

offspring_alpha = np.zeros((pop_size))	#----
offspring_beta = np.zeros((pop_size))	# Offspring
offspring_gamma = np.zeros((pop_size))	#
offspring_phi = np.zeros((pop_size))	#----

bin_alpha = np.zeros((pop_size - t), np.dtype('a6'))	#----
bin_beta = np.zeros((pop_size - t), np.dtype('a6'))		# Binary equivalents
bin_gamma = np.zeros((pop_size - t), np.dtype('a6'))	#
bin_phi = np.zeros((pop_size - t), np.dtype('a6'))		#----

pop_best = np.zeros((noi))
pop_avg = np.zeros((noi))
#---------------------------

#------------------
#	INITIALIZATION
#------------------
""" Randomly initialize the the four unknowns: the three euler angles 
	and the stress ellipsoid ratio """
alpha = np.random.randint(0, 64, size = pop_size)	#----
beta = np.random.randint(0, 64, size = pop_size)	# Euler angles alpha, beta and gamma
gamma = np.random.randint(0, 64, size = pop_size)	#-----
phi = np.random.randint(0, 64, size = pop_size)		# Stress ellipsoid ratio phi 

#-------------------------
#	ITERATIONS- MAIN CODE
#-------------------------
for n_i in range(0,noi):
	
	#	ENCODING
	""" Mapping the randomly initialized integers to real values in radians. Phi is just a ratio
		between 0 and 1."""
	alpha_en = 2 * np.pi * alpha / 63.00 
	beta_en = 2 * np.pi * beta / 63.00
	gamma_en = 2 * np.pi * gamma / 63.00
	phi_en = phi / 63.00
	
	#	GENERATING STRESS TENSORS
	""" Using the above aloha, beta, gamma and phi; we obtain the stress tensor in geographic 
		coordinate system."""
	for i in range(0, pop_size):
		reduced_tensor[i,:,:] = [[1, 0, 0],
								[0, phi_en[i], 0],
								[0, 0, 0]]

		t_matrix[i,:,:] = euler_rotation(euler_x(alpha_en[i]), euler_y(beta_en[i]), euler_z(gamma_en[i]))

		stress_tensor[i,:,:] = (np.transpose(t_matrix[i,:,:]).dot(reduced_tensor[i,:,:])).dot(t_matrix[i,:,:])

	#	FITNESS EVALUATION
	fitness = fitnessFunc(stress_tensor, pop_size, n_i)
	pop_best[n_i] = np.nanmin(fitness)	# Best fitness, i.e. minimum in the population
	pop_avg[n_i] = np.average(fitness)	# Average fitness of the population

	#	ELITISM
	""" The top t percent individuals do not undergo reproduction. The ensures that we 
		do not lose our best fit members in the next generation."""
	elite_index = np.argsort(fitness, axis = 0)
	for i in range(0, t):
		offspring_alpha[i] = alpha[elite_index[i]]
		offspring_beta[i] = beta[elite_index[i]]
		offspring_gamma[i] = gamma[elite_index[i]]
		offspring_phi[i] = phi[elite_index[i]]


	#	TOURNAMENT SELECTION
	""" Selection criteria for mating pool."""
	for i in range(0,pop_size - t):
		t1 = np.random.randint(t,pop_size)
		t2 = np.random.randint(t,pop_size)
		if fitness[t1] < fitness[t2]:
			m_alpha[i] = alpha[t1]
			m_beta[i] = beta[t1]
			m_gamma[i] = gamma[t1]
			m_phi[i] = phi[t1]
		else:
			m_alpha[i] = alpha[t2]
			m_beta[i] = beta[t2]
			m_gamma[i] = gamma[t2]
			m_phi[i] = phi[t2]

	#	CONVERTING DECIMAL TO BINARY
	for i in range(0, pop_size - t):
		bin_alpha[i] = dec2bin(int(m_alpha[i]), bits)
		bin_beta[i] = dec2bin(int(m_beta[i]), bits)
		bin_gamma[i] = dec2bin(int(m_gamma[i]), bits)
		bin_phi[i] = dec2bin(int(m_phi[i]), bits)


	#	CROSSOVER
	xover_alpha = crossoverFunc(bin_alpha, pop_size - t, bits)
	xover_beta = crossoverFunc(bin_beta, pop_size - t, bits)
	xover_gamma = crossoverFunc(bin_gamma, pop_size - t, bits)
	xover_phi = crossoverFunc(bin_phi, pop_size - t, bits)

	#	MUTATION
	mut_alpha = mutationFunc(xover_alpha, pop_size - t, bits)
	mut_beta = mutationFunc(xover_beta, pop_size - t, bits)
	mut_gamma = mutationFunc(xover_gamma, pop_size - t, bits)
	mut_phi = mutationFunc(xover_phi, pop_size - t, bits)


	for i in range(t, pop_size):
		offspring_alpha[i] = int(mut_alpha[i - t], 2)
		offspring_beta[i] = int(mut_beta[i - t], 2)
		offspring_gamma[i] = int(mut_gamma[i - t], 2)
		offspring_phi[i] = int(mut_phi[i - t], 2)

	#	UPDATING THE INITIAL POPULATION
	alpha = offspring_alpha
	beta = offspring_beta
	gamma = offspring_gamma
	phi = offspring_phi
#-------------------------------------------

#	THE BEST FIT INDIVIDUAL IN THE POPULATION
count = 0
for i in range(0,pop_size):
	if fitness[count] == np.nanmin(fitness):
		count = i
		break

for i in range(0,3):
	a1 = az_pl(t_matrix[count,i])
	azimuth = a1[0] if a1[0]>0 else (a1[0] + 360)
	plunge = a1[1]
	j = 3 - i 
	print "Sigma %d orientation is (Azimuth/Plunge):%d/%d " %(j,azimuth,plunge)

print "The final stress ellipsoid ratio is \n1. phi = %r" %(1 - phi_en[count])

print "Minimum Fitness" 
print np.nanmin(fitness)

print "Stress Tensor: " 
print stress_tensor[count,:,:]

x = np.arange(1,noi+1)
plt.plot(x, pop_avg, 'b-', x, pop_best, 'r-')
plt.xlabel('Iterations')
plt.ylabel('Misfit')
plt.show()	