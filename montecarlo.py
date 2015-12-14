
from scipy.stats import norm
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import math, random

class montecarlo(object):
    def __init__(self):
        self.Cover_Mean = 2.5 # cover depth in inches
        self.Cover_Stdev = .2  # cover standard deviation in inches
        self.Diffusivity_Mean = 0.005 # inches^2 / year
        self.Diffusivity_Stdev = 0.002  # inches^2 / year
        self.ChlorideThreshold_Mean = 2.6  # lbs / yd^3
        self.ChlorideThreshold_Stdev = 1.3  # lbs / yd^3
        self.ChlorideSurfaceConcentration_Mean = 25  # lbs / yd^3
        self.ChlorideSurfaceConcentration_Stdev = 10  # lbs / yd^3
        self.PropegationTime_Mean = 24  # crack propegation time years
        self.PropegationTime_Stdev = 12  # years
        self.Nitrite_addition = 0 #enter 1 for nitrite, 0 for no nitrite
        self.Nitrite_concentration = 10 #nitrite addition in pcy
        self.native_chloride = 0.1 #native chloride content in concrete, pcy
        self.crack_fraction = .01 #fraction of elements containing cracks
        self.crack_diffusivity_ratio = 10 #Diffusivity of crack containing elements relative to uncracked
        self.crack_diffusivity_stdev = 5
        self.crack_array = []
        self.jacket_cost = 10000 #installed Jacket cost
        self.jacket_replacement_time = 25 #jacket durability in years
        self.discount_rate = .05 #interest - inflation
        self.spalls_per_jacket = 3 #number of spalls showing before a jacket is required
        self.steps = 220000 #total number of elements in simulation
        self.derating_factor = [.4,.7] #derating factor i based on geometry, if you don't know what this is make it 1, allows for multiple deratings
        self.fraction_derated = [.4,.2] #fraction of elements with geometric derating i, allows for multiple, index is the same as derating factor's
        self.derating_array = []
        self.elements = 33 #elements per column
        self.time = 100 #time period ofver which the simulation is run
        self.sigma_range = 100#Truncation sigma range
        self.truncl = [1.1,.002,1,10,10,5]#truncation lower bounds for the 5 input parameters, the last number is for the crack diffusivity ratio
        self.trunch = [10,.1,20,60,60,20] #truncation upper bounds
        self.diffusion_type = 0 #0 for entered, 1 for rapid permeability, 2 for prediction based on concrete type
        self.rapid_permeability = 300 #output of the rapid permeability test in coulombs
        self.water_ratio = .3 #water to cement ratio
        self.cement_factor = 600 # lbs of cement in one cubic yard of concrete
        self.pozzolan_addition = 1 # are pozzolans present? enter 1 if yes, 0 if no
        self.element_calculation = 0 #how the program calculates the size of the simulation, 0 for manual entry, 1 for bridge piles
                #for a bridge deck or other flat section like a hotel front enter 2
        self.element_size = 3 #size of each discretized element in square feet, if no data assume 3
        self.column_perimeter = 10 #perimeter of each column in square feet, often ~10
        self.corrosion_zone = 10 #feet of area within corrosion zone
        self.number_of_piles = 2000 # often something like bridge length in feet divided by 15
        self.flat_footage = 200000 #total square footage of affected zone
        self.repair_size = 3 #size of section repaired, assume same as element size if repairs are simple patching
        self.rebar_type = 0 #0 for manual entry, 1 for plain steel, 2 for CRR, 3 for stainless
        self.derating_total = [0.0]
        self.monte_in = []
        self.monte_out = []
        self.monte_bincount = []
        self.cost = []
        self.jackets = []
        self.jacket_sum = [0.0]
        self.jacket_current = [0.0]
        self.cost_sum = [0.0]
        self.x = [0]
        self.b = 0


    def diffusion_tree(self):
        # This section calculates Diffusivity based on common user inputs
        if self.diffusion_type == 1:
            self.Diffusivity_Mean = .0005 * rapid_permeability**.78
            self.Diffusivity_Stdev = Diffusivity_Mean/2
        elif self.diffusion_type == 2:
            if self.pozzolan_addition == 1: F = 1
            else: F = 3
            self.Diffusivity_Mean = .01*(1+(self.water_ratio-.32)/.09)*(1+(752-self.cement_factor)/94)*F
            self.Diffusivity_Stdev = self.Diffusivity_Mean/2

    def structure_geometry(self):
        if self.element_calculation ==1:
            self.elements = int(self.column_perimeter * self.corrosion_zone/self.element_size)
            self.steps = self.elements * self.number_of_piles
        elif self.element_calculation ==2:
            self.steps = int(self.flat_footage/self.element_size)
            self.elements = self.repair_size / self.element_size

    def rebar_calc(self):
        if self.rebar_type == 1:
            self.ChlorideThreshold_Mean = .004 * self.cement_factor
        elif self.rebar_type == 2:
            self.ChlorideThreshold_Mean = .004 * self.cement_factor * 4
        elif self.rebar_type == 3:
            self.ChlorideThreshold_Mean = .004 * self.cement_factor * 10
        if self.rebar_type != 0:
            self.ChlorideSurfaceConcentration_Stdev = self.ChlorideThreshold_Mean/2

    def derating_calc(self):
        self.fraction_derated.append(1-sum(self.fraction_derated))
        self.derating_factor.append(1)
        for i in range(len(self.derating_factor)):
            self.derating_total.append(self.derating_total[i] + self.fraction_derated[i])
        self.derating_array = []
        for i in range(self.steps):
            for j in range(len(self.derating_factor)):
                if (i >=self.derating_total[j] * self.steps and i < self.derating_total[j+1] * self.steps):
                    self.derating_array.append(self.derating_factor[j])

    def crack_calc(self):
        for i in range (self.steps):
            t = 0
            while t == 0:
                proxy  = random.gauss(0,1)
                if self.truncl[5] < self.crack_diffusivity_ratio + proxy * self.crack_diffusivity_stdev < self.trunch[5]:
                    self.crack_array.append(self.crack_diffusivity_ratio + self.crack_diffusivity_stdev * proxy)
                    t = 1

    def monte_zeroes(self):
        for i in range(self.steps):
            self.monte_out.append(int(0))
        for i in range(5):
            t = []
            for j in range(self.steps):
                t.append(0.0)
            self.monte_in.append(t)

    def parameter_fill(self):
        self.averagevalues = [self.Cover_Mean, self.Diffusivity_Mean, self.ChlorideThreshold_Mean-self.native_chloride, self.ChlorideSurfaceConcentration_Mean,self.PropegationTime_Mean]
        self.StandardDeviations = [self.Cover_Stdev, self.Diffusivity_Stdev, self.ChlorideThreshold_Stdev, self.ChlorideSurfaceConcentration_Stdev,self.PropegationTime_Stdev]

    def monte_fill(self):
        if self.Nitrite_addition == 0:  #funct calculates the time to propegation for a given set of inputs
            def funct(a,b,c,d,e,f):
                if c/d < 1:
                    return np.square(a) / (4 * b * np.square(special.erfinv(1 - c / d))) + e #standard calculation of crack initiation time
                return 1000
        else:
            def funct(a,b,c,d,e,f):
                if ((f * (d-c)/(f+d)+c)/d) < 1:
                    return np.square(a) / (4 * b * np.square(special.erfinv(1 - (f * (d-c)/(f+d)+c)/d))) + e #crack initiation using Ct depentend on Nitrite addition and Cs
                return 1000
        for i in range(5):
            for j in range(self.steps):
                t = 0
                while t==0:
                    self.monte_in[i][j] = random.gauss(0,1)
                    a = 1
                    if i ==1:
                        a /= self.derating_array[j]
                    if random.random() < self.crack_fraction and i == 1: #accounting for cracks
                        a = self.crack_array[i]
                    if self.trunch[i] > self.monte_in[i][j] * self.StandardDeviations[i] + self.averagevalues[i] > self.truncl[i]:
                        self.monte_in[i][j] = self.monte_in[i][j] * self.StandardDeviations[i] + (a * self.averagevalues[i])
                        t=1
        for i in range(self.steps):
            self.monte_out[i]=funct(self.monte_in[0][i],self.monte_in[1][i],self.monte_in[2][i],self.monte_in[3][i],self.monte_in[4][i],self.Nitrite_concentration)
            self.monte_out[i]=math.floor(self.monte_out[i])

    def histogram(self):
        for i in range(int(self.time)):
            self.monte_bincount.append(0.0)
        for i in range(len(self.monte_out)):
            if len(self.monte_bincount) > self.monte_out[i] > 0:
                self.b = int(self.monte_out[i])
                self.monte_bincount[self.b] += 1

    def outputs(self):
        for i in range(0,len(self.monte_bincount)):
            if i > 0:
                self.monte_bincount[i] = self.monte_bincount[i]*((self.steps - self.elements*self.jacket_current[i-1])/self.steps)
            if i <= self.jacket_replacement_time - 1:
                self.jackets.append(self.monte_bincount[i]/self.spalls_per_jacket)
            else:
                self.jackets.append(self.monte_bincount[i]/self.spalls_per_jacket + self.jackets[i-self.jacket_replacement_time])
            if i > 0:
                self.jacket_sum.append(self.jackets[i] + self.jacket_sum[i-1])
                self.jacket_current.append(self.monte_bincount[i]/self.spalls_per_jacket + self.jacket_current[i-1])
            self.cost.append(self.jackets[i]*self.jacket_cost*((1-self.discount_rate)**i))
            if i>0:
                self.cost_sum.append(self.cost_sum[i-1]+self.cost[i])
    def normalization(self):
        self.monte_bincount[0] = self.monte_bincount[0]/self.steps
        for i in range(1,len(self.monte_bincount)):
            self.monte_bincount[i] = self.monte_bincount[i]/self.steps+self.monte_bincount [i-1]
            self.x.append(i)

    def plot_damage(self):
        plt.subplot(2,2,1)
        plt.plot(self.x,self.monte_bincount,'r--')
        plt.xlabel('Time (years)')
        plt.ylabel('fraction of elements showing spalls')
        plt.subplot(2,2,2)
        plt.plot(self.x,self.cost_sum,'r--')
        plt.xlabel('Time (years)')
        plt.ylabel('present cost of repair')
        plt.subplot(2,2,3)
        plt.plot(self.x,self.jacket_current,'r--')
        plt.xlabel('Time (years)')
        plt.ylabel('number of jackets currently on')
        plt.subplot(2,2,4)
        plt.plot(self.x,self.jacket_sum,'r--')
        plt.xlabel('Time (years)')
        plt.ylabel('total jackets installed')
        plt.show()

