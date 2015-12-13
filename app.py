# -*- coding: cp1252 -*-
from scipy.stats import norm
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import math, random

Cover_Mean = 2.5  # cover depth in inches
Cover_Stdev = .2  # cover standard deviation in inches
Diffusivity_Mean = 0.005 # inches^2 / year
Diffusivity_Stdev = 0.002  # inches^2 / year
ChlorideThreshold_Mean = 2.6  # lbs / yd^3
ChlorideThreshold_Stdev = 1.3  # lbs / yd^3
ChlorideSurfaceConcentration_Mean = 25  # lbs / yd^3
ChlorideSurfaceConcentration_Stdev = 10  # lbs / yd^3
PropegationTime_Mean = 24  # crack propegation time years
PropegationTime_Stdev = 12  # years
Nitrite_addition = 0 #enter 1 for nitrite, 0 for no nitrite
Nitrite_concentration = 10 #nitrite addition in pcy
native_chloride = 0.1 #native chloride content in concrete, pcy
crack_fraction = .01 #fraction of elements containing cracks
crack_diffusivity_ratio = 10 #Diffusivity of crack containing elements relative to uncracked
crack_diffusivity_stdev = 5
jacket_cost = 10000 #installed Jacket cost
jacket_replacement_time = 25 #jacket durability in years
discount_rate = .05 #interest - inflation
spalls_per_jacket = 3 #number of spalls showing before a jacket is required
steps = 220000 #total number of elements in simulation
derating_factor = [.4,.7] #derating factor i based on geometry, if you don't know what this is make it 1, allows for multiple deratings
fraction_derated = [.4,.2] #fraction of elements with geometric derating i, allows for multiple, index is the same as derating factor's
elements = 33 #elements per column
time = 100 #time period ofver which the simulation is run
sigma_range = 100#Truncation sigma range
truncl = [1.1,.002,1,10,10,5]#truncation lower bounds for the 5 input parameters, the last number is for the crack diffusivity ratio
trunch = [10,.1,20,60,60,20] #truncation upper bounds
# This section calculates Diffusivity based on common user inputs
diffusion_type = 0 #0 for entered, 1 for rapid permeability, 2 for prediction based on concrete type
rapid_permeability = 300 #output of the rapid permeability test in coulombs
water_ratio = .3 #water to cement ratio
cement_factor = 600 # lbs of cement in one cubic yard of concrete
pozzolan_addition = 1 # are pozzolans present? enter 1 if yes, 0 if no
if diffusion_type == 1:
    Diffusivity_Mean = .0005 * rapid_permeability**.78
    Diffusivity_Stdev = Diffusivity_Mean/2
elif diffusion_type == 2:
    if pozzolan_addition == 1: F = 1
    else: F = 3
    Diffusivity_Mean = .01*(1+(water_ratio-.32)/.09)*(1+(752-cement_factor)/94)*F
    Diffusivity_Stdev = Diffusivity_Mean/2
# end of section



#this section calculates the number of elements from common user iputs
element_calculation = 0 #how the program calculates the size of the simulation, 0 for manual entry, 1 for bridge piles
#for a bridge deck or other flat section like a hotel front enter 2
element_size = 3 #size of each discretized element in square feet, if no data assume 3
column_perimeter = 10 #perimeter of each column in square feet, often ~10
corrosion_zone = 10 #feet of area within corrosion zone
number_of_piles = 2000 # often something like bridge length in feet divided by 15
flat_footage = 200000 #total square footage of affected zone
repair_size = 3 #size of section repaired, assume same as element size if repairs are simple patching
if element_calculation ==1:
    elements = int(column_perimeter *corrosion_zone/element_size)
    steps = elements * number_of_piles
elif element_calculation ==2:
    steps = int(flat_footage/element_size)
    elements = repair_size / element_size
#end of section

#this section calculates the chloride threshold from common user inputs
rebar_type = 0 #0 for manual entry, 1 for plain steel, 2 for CRR, 3 for stainless
if rebar_type == 1:
    ChlorideThreshold_Mean = .004 * cement_factor
elif rebar_type == 2:
    ChlorideThreshold_Mean = .004 * cement_factor * 4
elif rebar_type == 3:
    ChlorideThreshold_Mean = .004 * cement_factor * 10
if rebar_type != 0:
    ChlorideSurfaceConcentration_Stdev = ChlorideThreshold_Mean/2
#end section

averagevalues = [Cover_Mean, Diffusivity_Mean, ChlorideThreshold_Mean-native_chloride, ChlorideSurfaceConcentration_Mean,
                 PropegationTime_Mean]
StandardDeviations = [Cover_Stdev, Diffusivity_Stdev, ChlorideThreshold_Stdev, ChlorideSurfaceConcentration_Stdev,
                      PropegationTime_Stdev]
# Assembles the above properties into tables of means and standard deviations

#Builds a table of deratings
derating_total = [0.0]
fraction_derated.append(1-sum(fraction_derated))
derating_factor.append(1)
for i in range(len(derating_factor)):
    derating_total.append(derating_total[i] + fraction_derated[i])
derating_array = []
for i in range(steps):
    for j in range(len(derating_factor)):
        if (i >=derating_total[j] * steps and i < derating_total[j+1] * steps):
            derating_array.append(derating_factor[j])
#End of section

#crack diffusivity variation
crack_array = []
for i in range (steps):
    t = 0
    while t == 0:
        proxy  = random.gauss(0,1)
        if truncl[5] < crack_diffusivity_ratio + proxy * crack_diffusivity_stdev < trunch[5]:
            crack_array.append(crack_diffusivity_ratio + crack_diffusivity_stdev * proxy)
            t = 1
#end of section
print crack_array[1]

if Nitrite_addition == 0:  #funct calculates the time to propegation for a given set of inputs
    def funct(a,b,c,d,e,f):
        if c/d < 1: 
            return np.square(a) / (4 * b * np.square(special.erfinv(1 - c / d))) + e #standard calculation of crack initiation time
        return 1000
else:
    def funct(a,b,c,d,e,f):
        return np.square(a) / (4 * b * np.square(special.erfinv(1 - (f * (d-c)/(f+d)+c)/d))) + e #crack initiation using Ct depentend on Nitrite addition and Cs
#here we are just filling in the matrices that we will use for the monte carlos simulation
monte_in = []
monte_out = []
monte_bincount = []

for i in range(steps):
    monte_out.append(int(0))

for i in range(5):
    t = []
    for j in range(steps):
        t.append(0.0)
    monte_in.append(t)
#the next sextion creates normally distributed variables with mean 0 and std 1 and multiplies them with the standard deviations and
#adds the mean, this gives a random value probibalistically distributed around the input mean
#there are two diffusivity corrections, basically the mean is multiplied by the crack diffusivity ratio or derating factor
#with some random probability equal to the input probability that a given element will have that characteristic
#monte_in is now populated with variables for each of the 5 input parameters, i references the input parameter and
#j references the element
for i in range(5):
    for j in range(steps):
        t = 0
        while t==0:
            monte_in[i][j] = random.gauss(0,1)
            a = 1
            if i ==1:
                a /= derating_array[j]
            if random.random() < crack_fraction and i == 1: #accounting for cracks
                a = crack_array[i]
            if  trunch[i] > monte_in[i][j] * StandardDeviations[i] + averagevalues[i] > truncl[i]:
                monte_in[i][j] = monte_in[i][j] * StandardDeviations[i] + (a * averagevalues[i])
                t=1
#the next section calculates a time to corrosion in years for each element in the monte_in array
#Monte_out is now populated with a time to corrosion for each element i, math.floor converts each element in
#monte_out to the integer below, since the outputs of funct are floats and we are going to be calculating corrosion times
#with a granularity of one year
for i in range(steps):
    monte_out[i]=funct(monte_in[0][i],monte_in[1][i],monte_in[2][i],monte_in[3][i],monte_in[4][i],Nitrite_concentration)
    monte_out[i] = math.floor(monte_out[i])
for i in range(int(time)):
    monte_bincount.append(0.0)
#the next section populates an array monte_mincount with the number failures that happened at any given year
#monte_bincount is basically a histogram of the times to corrosion, where each element i in the array references
#year i in the simulation
for i in range(len(monte_out)):
    if int(len(monte_bincount) > monte_out[i] > 0):
        a = int(monte_out[i])
        monte_bincount[a] += 1
cost = []
jackets = []
jacket_sum = [0.0]
jacket_current = [0.0]
cost_sum = [0.0]
#the next section calculates the number of elements damaged based on the number of jackets and the cathdic protection effect
#the number of jackets needed at any given year
#the number of jackets that have been put on total
#the number of jackets that are currently on
#the net present cost of repair for any given year
#the total net present cost of repair up until any given year
#each time a repair is given to a column the rest of the column is taken out of the simulation statistically
# such that the damaged fraction of the concrete converges to (damaged area per repair / size of each column)
for i in range(0,len(monte_bincount)):
    if i > 0:
        monte_bincount[i] = monte_bincount[i]*((steps - elements*jacket_current[i-1])/steps)
    if i <= jacket_replacement_time - 1:
        jackets.append(monte_bincount[i]/spalls_per_jacket)
    else:
        jackets.append(monte_bincount[i]/spalls_per_jacket + jackets[i-jacket_replacement_time])
    if i > 0: 
        jacket_sum.append(jackets[i] + jacket_sum[i-1])
        jacket_current.append(monte_bincount[i]/spalls_per_jacket + jacket_current[i-1])
    cost.append(jackets[i]*jacket_cost*((1-discount_rate)**i))
    if i>0:
        cost_sum.append(cost_sum[i-1]+cost[i])
#this just prints the numerical value of the present cost of maintenance at the service life of the structure
print "net present cost of total service life repair:"
print int(cost_sum[time-1])
x = [0]
monte_bincount[0]=monte_bincount[0]/steps
for i in range(1,len(monte_bincount)):
    monte_bincount[i] = monte_bincount[i]/steps+monte_bincount [i-1]
    x.append(i)
#here we print out our four plots
plt.subplot(2,2,1)
plt.plot(x,monte_bincount,'r--')
plt.xlabel('Time (years)')
plt.ylabel('fraction of elements showing spalls')

plt.subplot(2,2,2)
plt.plot(x,cost_sum,'r--')
plt.xlabel('Time (years)')
plt.ylabel('present cost of repair')

plt.subplot(2,2,3)
plt.plot(x,jacket_current,'r--')
plt.xlabel('Time (years)')
plt.ylabel('number of jackets currently on')

plt.subplot(2,2,4)
plt.plot(x,jacket_sum,'r--')
plt.xlabel('Time (years)')
plt.ylabel('total jackets installed')


plt.show()#outputs a graph of the damage function versus time
