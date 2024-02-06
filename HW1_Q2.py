import numpy as np
import math
from Bubble_dew import bubble_point_temperature,dew_point_temperature# Module defined by me to do Dew and Bubble point calculation based on Antonie eq and Raoults law
from Flash_code import flash_case1
names = ["methane", "ethylene", "propylene", "diethyl-ether", "ethanol", "isopropanol", "water"]
coefficients = list([[6.61184, 389.93, 266.0],\
                      [6.74756, 585.,255.],\
                      [6.81960, 785., 247.],\
                      [7.4021,1391.4, 273.16],\
                      [8.04494,1554.3,222.65],\
                      [6.66040, 813.055, 132.96],\
                      [8.10765,1750.286,235.]])
# unpacking the coefficients list A,B,C to be used in future calculations.
A= [sub_array[0] for sub_array in coefficients]
B= [sub_array[1] for sub_array in coefficients]
C= [sub_array[2] for sub_array in coefficients]
fk=np.array([200.0,1198.77,266.71,2.421,90.79,1.8802,680.72]) # feed is in vapour phase as per Table 3.1
#fk =fk *3.6 #----> Convert to kgmol/hr

yi=np.array(fk/np.sum(fk))
Temp_feed=590.0-273.15 # input stream temperature for flash , Kelvin--> celcius
Press_operating= 68.5*750.062 # in bar to mmHg

# we have to note the calculation for antonie eq will be done in mmHg and Celcius

'''Answer'''
print('\n Answer 1\n')
dew_point_temperature(P=Press_operating,A=A,B=B,C=C,y=yi)
bubble_point_temperature(P=Press_operating,A=A,B=B,C=C,x=yi)
print('\n Answer 2')
flash_case1(A=A,B=B,C=C,fk=fk,x=False)# pressure specified so x = False









