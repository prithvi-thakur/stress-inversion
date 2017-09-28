# stress-inversion
Python 2.7 programs for forward modelling and inversion of stress from fault-slip observations using genetic algorithm.

There are four files-

SyntheticData.py: creates a set of 20 ideal synthetic fault-slip data sets. You can manipulate the values of sigma1, sigma3 and the stress ratio (phi) values to get the desired stress regime. Running this program creates an excel sheet with its columns as dip direction, dip, trend of lineation and the plunge of lineation. Note that the plunge will be negative for reverse faults and positive for normal faults. The reference frame is North, East and Up as x, y and z respectively.
  
demo.xls: This is the excel sheet of synthetic data sets. Just run the synthetic data with the desired stress tensor as the input to change the contents of this file.
  
GAstress.py: This is the stress inversion program for homogeneous fault-slip data. It is based on the genetic algorithms to perform inversion of fault-slip data to get the reduced stress tensor. For details refer: http://www.sciencedirect.com/science/article/pii/S0191814116302036

hga.py: The stress inversion program for heterogeneous fault-slip data. The algorithm is described in detail in my paper _____. I'll add it here once it is accepted.
  
  
SPECIFIC INSTRUCTIONS:
------------------
Input Data format:
------------------
The input data is an excel sheet with four columns: 

dip direction, dip angle of fault plane, azimuth, and plunge of lineation.

We follow a left handed North East Up coordinate system where 
Dip direction is 90 degrees clockwise from the strike.



--------------
Instructions:
--------------

1. The excel data file should be in the same directory as the python script.

2. Python must be installed in the system. (You can use IDEs like pycharm or directly install it from official python website)

3. Open the python file HGA.py or GAstress.py

4. Go to line 82(GAstress.py) or line 63(hga.py), and within the paranthesis change the name of the input data file as desired along with the file extension. E.g: "demo.xlsx", "data2.xls", etc.

5. <b>For running hga.py:</b> go to line 207, and change the value of M to the number of expected stress states. The default value is 4. It is advisable to run the program once with the value of M (step 5 above) more than the expected stress states and a second time with the value of M equal to the expected stress states.

6. Run the script. Wait for a few minutes before result pops up.







--
For any queries, email prithvi.thakur93@gmail.com
