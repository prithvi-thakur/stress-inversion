#	-------------------------------------------------------------
#	-------------------------------------------------------------
# 	 	THE GENETIC ALGORITHM METHOD FOR STRESS INVERSION
#		Python Script to invert homogeneous fault-slip data
#		Version: 1.0
#		Author: Prithvi Thakur
#		Last Modified: 14 November 2016, 7:20 p.m. 
#	-------------------------------------------------------------
#   -------------------------------------------------------------	


#	Copyright 2016 Prithvi Thakur 
#		
#	This program is free software: you can redistribute it and/or modify
#  	it under the terms of the GNU General Public License as published by
#  	the Free Software Foundation, either version 3 of the License, or
#  	(at your option) any later version.

# 	This program is distributed in the hope that it will be useful,
# 	but WITHOUT ANY WARRANTY; without even the implied warranty of
#  	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  	GNU General Public License for more details.

#  	You should have received a copy of the GNU General Public License
#  	along with this program.  If not, see <http://www.gnu.org/licenses/>.


#	INSTRUCTIONS:

#	Line 82 takes the input data in excel file format with four columns: 
#	dip direction, dip angle, trend and plunge of lineations. 
#	The input data file should be in the same directory as this
#	program. Just add the file name of the data you want to invert, and
#	run it in python.

#	You can play with the variables in lines 219-223 to tweak the algorithm
#	and possibly obtain better inversion results.


#	Import python packages
import numpy as np 
import random
import math
import xlrd
import matplotlib.pyplot as plt

random.seed()	#seed the random number - system time


#------------------------------------
#	ROTATION MATRIX	(check wikipedia)
#------------------------------------
euler_z = lambda z: np.array([[np.cos(z), -np.sin(z), 0], [np.sin(z), np.cos(z), 0] ,[0, 0, 1]])
euler_y = lambda y: np.array([[np.cos(y), 0, np.sin(y)], [0, 1, 0], [-np.sin(y), 0, np.cos(y)]])
euler_x = lambda x: np.array([[1, 0, 0], [0, np.cos(x), -np.sin(x)], [0, np.sin(x), np.cos(x)]])

euler_rotation = lambda x, y, z: z.dot(y.dot(x))


#-------------------------------------------
#	DIRECTION COSINES FROM AZIMUTH/PLUNGE
#-------------------------------------------
dir_cosine = lambda az, pl: np.array([np.cos(pl)*np.sin(az), np.cos(pl)*np.cos(az), -np.sin(pl)])

#-------------------------------------------
#	AZIMUTH/PLUNGE FROM DIRECTION COSINE
#-------------------------------------------
az_pl = lambda vector: np.array([np.degrees(np.arctan2(vector[0], vector[1])), 
								np.degrees(np.arcsin(-vector[2]/np.linalg.norm(vector)))])

#--------------------------------------------------
#	DECIMAL TO BINARY FUNCTION, n = number of bits
#--------------------------------------------------
dec2bin = lambda x, n:format(x, 'b').zfill(n) 


#---------------------
#	FITNESS FUNCTION
#---------------------
def fitnessFunc(stress_tensor, M, pop_size, n_i):
	"This function computes the composite misfit for M stress states"
	input_data = xlrd.open_workbook('N2.xls')
	first_sheet = input_data.sheet_by_index(0)
	dip_d = np.radians(first_sheet.col_values(0))	#	Dip Direction
	dip_a = np.radians(first_sheet.col_values(1))	#	Dip Angle
	trend = np.radians(first_sheet.col_values(2))	#	Trend of Lineation
	plunge = np.radians(first_sheet.col_values(3))	#	Plunge of Lineation
	length = len(dip_d)

	#	PREALLOCATION
	normal_vector = np.zeros((length,3))	#	Normal Vector
	slip_vector = np.zeros((length,3))		#	Slip Vector
	fitness = np.zeros((pop_size,1))		#	Misfit value
	A_array = [np.zeros((length)) for x in range(M)]
	DMM = [np.zeros((length)) for x in range(M)]

	#	SUB FITNESS FUNCTION
	def subfitnessFunc(sigma, slip, normal):
		"This function computes the composite misfit for each phase, each datum"
		stress_vector = sigma.dot(normal)
		observed_shear_component = stress_vector.dot(slip)
		stress_vector_mag = np.linalg.norm(stress_vector, 2)
		normal_component = stress_vector.dot(normal)
		shear_component = stress_vector - (normal * normal_component)
		unit_shear_component = shear_component / np.linalg.norm(shear_component)

		new_misfit = np.absolute(np.arccos(unit_shear_component.dot(slip)))

		f2 = np.sin((new_misfit/2))
		f1 = np.absolute( np.absolute(observed_shear_component) - np.sqrt(np.square(stress_vector_mag) - np.square(normal_component)))

		lamb = 0.7
		f = lamb*f1 + (1-lamb)*f2

		return f

	#	COMPUTING NORMAL AND SHEAR VECTOR
	for i in range(0, length):	
		normal_vector[i,:] = dir_cosine(dip_d[i], dip_a[i]-(np.pi/2))
		slip_vector[i,:] = dir_cosine(trend[i], plunge[i])

	#	COMPUTING FITNESS ARRAY
	for q in range(0,pop_size):
		sum_f = [0 for x in range(M)]
		f = [0 for x in range(M)]
		counter = [0 for x in range(M)]
		mean_f = [0 for x in range(M)]
		separated = np.zeros((length))

		for w in range(0,length):
			normal = normal_vector[w,:]
			slip = slip_vector[w,:]

			for i1 in range(M):
				f[i1] = subfitnessFunc(stress_tensor[i1][q,:,:], slip, normal)
				DMM[i1][w] = f[i1]

			for i1 in range(M):
				if f[i1] == min(f):
					sum_f[i1] = sum_f[i1] + f[i1]
					separated[w] = i1 + 1
					counter[i1]+= 1
					A_array[i1][w] = f[i1]

		for i1 in range(M):
			mean_f[i1] = 1 if (counter[i1] == 0) else (sum_f[i1]/counter[i1])

		sum_f = sum(mean_f)

		for i1 in range(M):
			A_ar = filter(lambda x: x!=0, A_array[i1])

		fitness[q] = sum_f #+ std_f

	if n_i == noi-1:
		#print "The number of data points in stress tensor: \n" 
		#for i1 in range(M):
		#	print counter[i1]
		#print separated
		print "Minimum fitness"
		print np.nanmin(fitness)


	return fitness
#--------------------	


#---------------------
#	CROSSVER FUNCTION
#---------------------
def crossoverFunc(parents, size, bits):
	"""This function performs single point crossover for input array of strings.
	The crossover probability is 60 percent"""
	children = np.zeros((size), np.dtype('a6'))

	for i in range(0, int(size/2)):
		x_site = np.random.randint(0, bits - 1)
		x1 = parents[i]
		x2 = parents[size - i - 1]
		if (np.random.randint(0, 100)) > 40 :	# Crossover Probability = 60 percent
			ch1 = x1[0:x_site] + x2[x_site:bits]
			ch2 = x2[0:x_site] + x1[x_site:bits]

			children[i] = ch1
			children[size - i - 1] = ch2
		
		else:
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

	if prob > 90:	# Mutation Probability = 90 percent
		s = np.random.randint(0, size)
		b = np.random.randint(0, bits)
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
pop_size = 800		# Population size
noi = 60			# Number of Iterations
bits = 6			# Number of bits of each parameter
M = 1				# Number of Phases in the Data
elite = 10			# Percentage of elite individuals
t = int(elite * 0.01 * pop_size)	# Fraction of elite individuals 
print "The Population size is %d and the number of iterations is %d" % (pop_size, noi)

reduced_tensor = [np.zeros((pop_size,3,3)) for x in range(M)]
stress_tensor = [np.zeros((pop_size,3,3)) for x in range(M)]
t_matrix = [np.zeros((pop_size,3,3)) for x in range(M)]

fitness = np.zeros((pop_size))

m_alpha = [np.zeros((pop_size-t)) for x in range(M)]
m_beta = [np.zeros((pop_size-t)) for x in range(M)]
m_gamma = [np.zeros((pop_size-t)) for x in range(M)]
m_phi = [np.zeros((pop_size-t)) for x in range(M)]

offspring_alpha = [np.zeros((pop_size)) for x in range(M)]
offspring_beta = [np.zeros((pop_size)) for x in range(M)]
offspring_gamma = [np.zeros((pop_size)) for x in range(M)]
offspring_phi = [np.zeros((pop_size)) for x in range(M)]


bin_alpha = [np.zeros((pop_size - t), np.dtype('a6')) for x in range(M)]
bin_beta = [np.zeros((pop_size - t), np.dtype('a6')) for x in range(M)]
bin_gamma = [np.zeros((pop_size - t), np.dtype('a6')) for x in range(M)]
bin_phi = [np.zeros((pop_size - t), np.dtype('a6')) for x in range(M)]

xover_alpha = [np.zeros((pop_size - t)) for x in range(M)]
xover_beta = [np.zeros((pop_size - t)) for x in range(M)]
xover_gamma = [np.zeros((pop_size - t)) for x in range(M)]
xover_phi = [np.zeros((pop_size - t)) for x in range(M)]

mut_alpha = [np.zeros((pop_size - t)) for x in range(M)]
mut_beta = [np.zeros((pop_size - t)) for x in range(M)]
mut_gamma = [np.zeros((pop_size - t)) for x in range(M)]
mut_phi = [np.zeros((pop_size - t)) for x in range(M)]

pop_best = np.zeros((noi))
pop_avg = np.zeros((noi))
#---------------------------


#------------------
#	INITIALIZATION
#------------------
alpha = [np.random.randint(0, 64, size = pop_size) for x in range(M)]
beta = [np.random.randint(0, 64, size = pop_size) for x in range(M)]
gamma = [np.random.randint(0, 64, size = pop_size) for x in range(M)]
phi = [np.random.randint(0, 64, size = pop_size) for x in range(M)]

#-------------------------
#	ITERATIONS- MAIN CODE
#-------------------------
for n_i in range(0,noi):
	
	#	ENCODING
	alpha_en = [(2*np.pi*x/63.00) for x in alpha]
	beta_en = [((np.pi/2)*x/63.00) for x in beta]
	gamma_en = [(2*np.pi*x/63.00) for x in gamma]
	phi_en = [(x/63.00) for x in phi]
	
	#	GENERATING STRESS TENSORS
	for i in range(M):
		for j in range(pop_size):
			reduced_tensor[i][j,:,:] = [[1, 0, 0],
								[0, phi_en[i][j], 0],
								[0, 0, 0]]

			t_matrix[i][j,:,:] = euler_rotation(euler_x(alpha_en[i][j]), euler_y(beta_en[i][j]), euler_z(gamma_en[i][j]))
			stress_tensor[i][j,:,:] = (t_matrix[i][j].dot(reduced_tensor[i][j])).dot(np.transpose(t_matrix[i][j]))

	#	FITNESS EVALUATION
	fitness = fitnessFunc(stress_tensor, M, pop_size, n_i)
	pop_best[n_i] = np.nanmin(fitness)
	pop_avg[n_i] = np.average(fitness)

	#	ELITISM
	elite_index = np.argsort(fitness, axis = 0)
	for i in range(M):
		for j in range(t):
			offspring_alpha[i][j] = alpha[i][elite_index[j]]
			offspring_beta[i][j] = beta[i][elite_index[j]]
			offspring_gamma[i][j] = gamma[i][elite_index[j]]
			offspring_phi[i][j] = phi[i][elite_index[j]]

	#	TOURNAMENT SELECTION
	for i in range(M):
		for j in range(0,pop_size - t):
			t1 = np.random.randint(t,pop_size)
			t2 = np.random.randint(t,pop_size)
			if fitness[t1] < fitness[t2]:
				m_alpha[i][j] = alpha[i][t1]
				m_beta[i][j] = beta[i][t1]
				m_gamma[i][j] = gamma[i][t1]
				m_phi[i][j] = phi[i][t1]

			else:
				m_alpha[i][j] = alpha[i][t2]
				m_beta[i][j] = beta[i][t2]
				m_gamma[i][j] = gamma[i][t2]
				m_phi[i][j] = phi[i][t2]

	#	CONVERTING DECIMAL TO BINARY
	for i in range(M):
		for j in range(0, pop_size - t):
			bin_alpha[i][j] = dec2bin(int(m_alpha[i][j]), bits)
			bin_beta[i][j] = dec2bin(int(m_beta[i][j]), bits)
			bin_gamma[i][j] = dec2bin(int(m_gamma[i][j]), bits)
			bin_phi[i][j] = dec2bin(int(m_phi[i][j]), bits)

	#	CROSSOVER
	for i in range(M):
		xover_alpha[i] = crossoverFunc(bin_alpha[i], pop_size - t, bits)
		xover_beta[i] = crossoverFunc(bin_beta[i], pop_size - t, bits)
		xover_gamma[i] = crossoverFunc(bin_gamma[i], pop_size - t, bits)
		xover_phi[i] = crossoverFunc(bin_phi[i], pop_size - t, bits)

	#	MUTATION
	for i in range(M):
		mut_alpha[i] = mutationFunc(xover_alpha[i], pop_size - t, bits)
		mut_beta[i] = mutationFunc(xover_beta[i], pop_size - t, bits)
		mut_gamma[i] = mutationFunc(xover_gamma[i], pop_size - t, bits)
		mut_phi[i] = mutationFunc(xover_phi[i], pop_size - t, bits)

	for i in range(M):
		for j in range(t, pop_size):
			offspring_alpha[i][j] = int(mut_alpha[i][j - t], 2)
			offspring_beta[i][j] = int(mut_beta[i][j - t], 2)
			offspring_gamma[i][j] = int(mut_gamma[i][j - t], 2)
			offspring_phi[i][j] = int(mut_phi[i][j - t], 2)

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

for i in range(M):
	for j in range(0,3):
		a1 = az_pl(t_matrix[i][count,j,:])
		azimuth = a1[0] if a1[0]>0 else (a1[0] + 360)
		if a1[1] < 0:
			azimuth = (azimuth -180) if (azimuth-180)>0 else (azimuth + 180)
			plunge = -a1[1]
		else:
			plunge = a1[1]
		print "Sigma %d of stress regime %d is:%d/%d " %(3-j, i+1, azimuth,plunge)

for i in range(M):
	print "Stress Ratio of stress regime %d is: %r" %(i+1, 1-phi_en[i][count])

x = np.arange(1,noi+1)
plt.plot(x, pop_avg, 'b-', x, pop_best, 'r-')
plt.xlabel('Iterations')
plt.ylabel('Misfit')
plt.show()	