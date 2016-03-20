# stress-inversion
Python program for forward modelling and inversion of stress from fault-slip observations using genetic algorithms

There are three files-

SyntheticData: creates a set of 20 ideal synthetic fault-slip data sets. You can manipulate the values of sigma1, sigma3
  and the stress ratio (phi) values to get the desired stress regime. Running this program creates an excel sheet with its
  columns as dip direction, dip, trend of lineation and the plunge of lineation. Note that the plunge will be negative for
  reverse faults and positive for normal faults. The reference frame is North, East and Up as x, y and z respectively.
  
Uselessdata: This is the excel sheet of synthetic data sets. Just run the synthetic data with the desired stress tensor as
  the input to change the contents of this file.
  
Inversion: This is the actual inversion program. It is based on the genetic algorithms to perform inversion of fault-slip
  data to get the reduced stress tensor. The input data has to be in excel file format in the same directory as this program.
  Just change the file name in line 55 to change the input data set. The generalized theory og this inversion can be found 
  in Angelier (1994) paper. I've submitted a paper for my algorithm, and will refer it here once it is accepted.
  
  
