#!/usr/bin/env python
# coding: utf-8

# In[ ]:


print("Bayode")


# We are currently in the first class for Computational Physics. Basically trying things out.         

# Exercise 2.1

# In[ ]:


import math
import numpy as np


# In[ ]:


h = float(input("Enter the height of the tower: "))
t = float(input("Enter the time interval: "))

s = 9.81*t**2/2

print("The height of the ball is",h-s,"meters")


# Enter the height of the tower: 100
# Enter the time interval: 5
# The height of the ball is -22.625 meters

# In[ ]:


h = float(input("Enter the height of the tower: "))
g = 9.81

s = 9.81*t**2/2

t = (2*(g**-1)*h)**(1/2)

print("The time is",t)


# # Exercise 2.2: Altitude of a satellite
# 
# The formula for calculating the altitude of a satellite in a circular orbit around the Earth is:
# 
# $$h = \left(\frac{G \cdot M \cdot T^2}{4 \cdot \pi^2}\right)^{\frac{1}{3}} - R$$
# 
# Where:
# h is the altitude of the satellite above the Earth's surface (in meters)
# G is the gravitational constant (approximately 6.67 x 10^-11 N*(m^2)/(kg^2))
# M is the mass of the Earth (approximately 5.97 x 10^24 kg)
# T is the orbital period of the satellite (in seconds)
# pi is approximately 3.14
# R is the radius of the Earth (approximately 6,371,000 meters)
# 
# This formula is known as "Kepler's Third Law" which states that the square of the orbital period of a planet is proportional to the cube of the semi-major axis of its orbit.

# In[ ]:


###Exercise 2:2 b

import math

G = 6.67 * (10 ** -11)
M = 5.97 * (10 ** 24)
R = 6371 * (10 ** 6)

T = float(input("Enter the desired orbital period in seconds: "))

h = ((G * M * T ** 2) / (4 * math.pi ** 2)) ** (1/3) - R

print("The altitude of the satellite in meters is: ", h)


# 

# In[9]:


###Exercise 2:2 c
    

import math

G = 6.67 * (10 ** -11)
M = 5.97 * (10 ** 24)
R = 6371 * (10 ** 3)

T = 86400
h = ((G * M * T ** 2) / (4 * math.pi ** 2)) ** (1/3) - R
print("The altitude of the satellite in a geosynchronous orbit in meters is: ", h)

T = 5400
h = ((G * M * T ** 2) / (4 * math.pi ** 2)) ** (1/3) - R
print("The altitude of the satellite in a 90 minutes orbit in meters is: ", h)

T = 2700
h = ((G * M * T ** 2) / (4 * math.pi ** 2)) ** (1/3) - R
print("The altitude of the satellite in a 45 minutes orbit in meters is: ", h)


# A geosynchronous orbit is an orbit around the Earth with an orbital period of exactly 24 hours (86400 seconds) and it is used for communication satellites.
# A 90 minutes orbit is a Low Earth Orbit (LEO) and it's used for remote sensing, and earth observation.
# A 45 minutes orbit is a very low Earth orbit, which is not a practical orbit for any known satellite and this is because it's too close to the Earth's surface, it is affected by drag and atmospheric forces, and it's also too close to the Van Allen radiation belts, making it not suitable for any kind of payload.
# Hence, it is concluded that a 45 minutes orbit is not practical for a satellite to be launched into.

# ### Exercise 2:2 d
# 
# A geosynchronous satellite is one that orbits the Earth once per sidereal day which is 23.93 hours, not 24 hours. This is because the Earth is not only rotating on its axis but also orbiting around the sun. The rotation of the Earth takes 24 hours to complete, but because the Earth is also orbiting the sun, it takes an additional 4 minutes (about 1/365th of a day) for the Earth to return to the same position in its orbit. Therefore, a geosynchronous satellite must orbit the Earth once every 23 hours and 56 minutes, or 86,164 seconds, in order to keep pace with the rotation of the Earth and remain above the same point on the surface.
# 
# The difference in altitude between a geosynchronous satellite orbiting once every 24 hours and a geosynchronous satellite orbiting once every 23.93 hours is very small.
# It can be calculated by using the formula provided in the previous answer but substituting 86400 seconds with 86164 seconds.
# It will make a difference of about 36 meters to the altitude of the satellite.
# This small difference is not significant in practical terms and geosynchronous orbit is still commonly defined as an orbit with an orbital period of 24 hours.

# ### Exercise 2.3: 
# 
# Write a program to perform the inverse operation to that of Example 2.2. That is, ask the user for the Cartesian coordinates x, y of a point in two-dimensional space, and calculate and print the corresponding polar coordinates, with the angle 9 given in degrees.

# In[16]:


from math import sin, cos, pi, acos
from numpy import rad2deg
x = float(input("Enter x: "))
y = float(input("Enter y: "))

r = (x**2) + (y**2)
theta = rad2deg(acos(x/r))
print("The value for r is: ", r)
print("The value of theta is: ", theta )


# ### Exercise 2.4:
# 
# A spaceship travels from Earth in a straight line at relativistic speed v to another planet x light years away. Write a program to ask the user for the value of x and the speed v as a fraction of the speed of light c, then print out the time in years that the spaceship takes to reach its destination (a) in the rest frame of an observer on Earth and (b) as perceived by a passenger on board the ship. Use your program to calculate the answers for a planet 10 light years away with v = 0.99c.

# In[ ]:


import math

c = 3*10**8

x = float(input("Enter the distance to the planet in light years: "))
v = float(input("Enter the speed of the spaceship as a fraction of the speed of light (e.g. 0.99 for 0.99c): "))

# time in years as perceived by an observer on Earth
t_earth = x / (math.sqrt(1 - (v**2/c**2)) * v)

# time in years as perceived by a passenger on board the ship
t_ship = x / (v * math.sqrt(1 - (v**2/c**2)))

print("Time in years as perceived by an observer on Earth: ", t_earth)
print("Time in years as perceived by a passenger on board the ship: ", t_ship)


# ###Exercise 2.5: Quantum potential step
# 
# 
# A well-known quantum mechanics problem involves a particle of mass m that encounters a one-dimensional potential step, like this: (Image is textbook)
# 
# 
# The particle with initial kinetic energy E and wavevector $$k_1 = \sqrt{\frac{2mE}{h}}$$ enters from the left and encounters a sudden jump in potential energy of height V at position x = 0. By solving the Schrodinger equation, one can show that when E > V the particle may either (a) pass the step, in which case it has a lower kinetic energy of E - V on the other side and a correspondingly smaller wavevector of $$k_2 = \sqrt{\frac{2m(E - V)}{h}}$$, or 
# 
# (b) it may be reflected, keeping all of its kinetic energy and an unchanged wavevector but moving in the opposite direction. The probabilities T and R for transmission and reflection are given by
# 
# $$T = \frac{4k_1k_2}{(k_1 + k_2)^2}$$ , $$R = \left(\frac{k_1 - k_2}{k_1 + k_2}\right)^2$$
# 
# Suppose we have a particle with mass equal to the electron mass m = 9.11 x 10^-31 kg and energy 10 eV encountering a potential step of height 9eV. Write a Python program to compute and print out the transmission and reflection probabilities using the formulas above.
# 
# 

# In[49]:


import math 
import scipy.constants

m = 9.11*(10**(-31))
E = 10 * (1.6)*10**(-19)
V = 9 * (1.6)*10**(-19)

hbar = (6.582119569*10**(-16)) * (1.6)*10**(-19)

k1 = math.sqrt(2*m*E)/hbar
k2 = math.sqrt(2*m*(E - V))/hbar


T = (4*(k1)*(k2))/((k1)+(k2))**2

R = (((k1)-(k2))/((k1)+(k2)))**2



print("The value of T is: ", T)

print("The value of R is: ", R)


# ### Exercise 2.6: Planetary orbits
# 
# The orbit in space of one body around another, such as a planet around the Sun, need not be circular. In general it takes the form of an ellipse, with the body sometimes closer in and sometimes further out. If you are given the distance ℓ1 of closest approach that a planet makes to the Sun, also called its perihelion, and its linear velocity v1 at perihelion, then any other property of the orbit can be calculated from these two as follows.
# 
# a. Kepler's second law tells us that the distance l2 and velocity v2 of the planet at its most distant point, or aphelion, satisfy $$l_2v_2 = l_1v_1$$• At the same time the total the total energy, kinetic plus gravitational, of a planet with velocity v and distance r from the Sun is given by
# $$E = \frac{1}{2}mv^2 - \frac{GmM}{r}$$
# where m is the planet’s mass, M = 1.9891 × 1030 kg is the mass of the Sun, and G = 6.6738 × 10^−11 m^3 kg^−1 s^−2 is Newton’s gravitational constant. Given that energy must be conserved, show that v2 is the smaller root of the quadratic equation.
# 
# $$v_2^2 - 2\frac{GM}{l_1}v_2\frac{v_1}{l_1} - (v_1^2 - 2\frac{GM}{l_1}) = 0$$
# 
# Once we have v2 we can calculate l2 using the relation $$l_2 = \frac{l_1 v_1}{v_2}$$.
# 
# b) Given the values of v1, l1, and l2, other parameters of the orbit are given by simple formulas can that be derived from Kepler’s laws and the fact that the orbit is an ellipse:
# 
# Semi-major axis: $$a = \frac{1}{2}(l_1 + l_2)$$
# 
# Semi-minor axis: $$b = \sqrt{l_1 l_2}$$
# 
# Orbital period: $$T = \frac{2\pi ab}{l_1 v_1}$$
# 
# Orbital eccentricity: $$e = \frac{l_2 - l_1}{l_2 + l_1}$$
# 
# 
# Write a program that asks the user to enter the distance to the Sun and velocity at perihelion, then calculates and prints the quantities l2, v2, T, and e.
# 
# c) Test your program by having it calculate the properties of the orbits of the Earth (for which l1 = 1.4710 × 10^11 m and v1 = 3.0287 × 10^4 m s−1) and Halley’s comet (l1 = 8.7830 × 10^10 m and v1 = 5.4529 × 10^4 m s−1). Among other things, you should find that the orbital period of the Earth is one year and that of Halley’s comet is about 76 years.

# In[ ]:


import math

l1 = float(input("Enter the distance to the Sun (in AU): "))
v1 = float(input("Enter the velocity at perihelion (in km/s): "))

l2 = l1 * v1 / v2
a = 0.5 * (l1 + l2)
b = math.sqrt(l1 * l2)
T = (2 * math.pi * a * b) / (l1 * v1)
e = (l2 - l1) / (l2 + l1)

print("l2: ", l2, "AU")
print("v2: ", v2, "km/s")
print("T: ", T, "years")
print("e: ", e)


# In[52]:


import math

M = 1.9891 * 10**30
G = 6.6738 * 10 **-11

# Earth's orbit
l1 = 1.4710*(10**11) # distance to the sun in m
v1 = 3.0287*(10**4) # velocity at perihelion in m/s

# Halley's comet
#l1 = 8.7830e10 # distance to the sun in m
#v1 = 5.4529e4 # velocity at perihelion in m/s

#Calculating the properties of the orbit

v2 = math.sqrt((2*G*M)/l1 + v1**2)
l2 = l1 * v1 / v2
a = 0.5 * (l1 + l2)
b = math.sqrt(l1 * l2)
T = (2 * math.pi * a * b) / (l1 * v1)
e = (l2 - l1) / (l2 + l1)

print("l2: ", l2/1000, "Km")
print("v2: ", v2, "m/s")
print("T: ", T/31536000, "years")
print("e: ", e)


# In[53]:


import math

M = 1.9891 * 10**30
G = 6.6738 * 10 **-11

# Earth's orbit
#l1 = 1.4710*(10**11) # distance to the sun in m
#v1 = 3.0287*(10**4) # velocity at perihelion in m/s

# Halley's comet
l1 = 8.7830e10 # distance to the sun in m
v1 = 5.4529e4 # velocity at perihelion in m/s

#Calculating the properties of the orbit

v2 = math.sqrt((2*G*M)/l1 + v1**2)
l2 = l1 * v1 / v2
a = 0.5 * (l1 + l2)
b = math.sqrt(l1 * l2)
T = (2 * math.pi * a * b) / (l1 * v1)
e = (l2 - l1) / (l2 + l1)

print("l2: ", l2/1000, "Km")
print("v2: ", v2, "m/s")
print("T: ", T/31536000, "years")
print("e: ", e)


# ### Exercise 2.7: Catalan numbers
# 
# The Catalan numbers Cn are a sequence of integers 1, 1, 2, 5, 14, 42, 132... that play an important role in quantum mechanics and the theory of disordered systems. (They were central to Eugene Wigner’s proof of the so-called semicircle law.) They are given by
# 
# $$C_0 = 1$$
# 
# $$C_{n+1} = \frac{(4n+2)}{(n+2)}C_n$$
# 
# Write a python program that prints in increasing order all Catalan numbers less than or equal to one billion.

# In[54]:


def catalan(n):
    if n == 0:
        return 1
    else:
        return int((4*n-2)/(n+1)*catalan(n-1))

n = 0
while catalan(n) <= 1000000000:
    print(catalan(n))
    n += 1


# This program uses a recursive function catalan(n) to calculate the nth Catalan number using the formula provided. The function has a base case of n = 0, which returns 1. For n > 0, the function returns ((4n-2)/(n+1))*catalan(n-1).
# 
# The program then enters a while loop, which continues as long as the value of catalan(n) is less than or equal to one billion. In each iteration of the loop, the current value of catalan(n) is printed and n is incremented.
# 
# Note that this program will print the catalan number and will not store them in a list or array. If you want to store them in an array you can use a list and append the catalan number to it.

# In[56]:


catalan_numbers = []
def catalan(n):
    if n == 0:
        return 1
    else:
        return int((4*n-2)/(n+1)*catalan(n-1))

n = 0
while catalan(n) <= 1000000000:
    catalan_numbers.append(catalan(n))
    n += 1

print(catalan_numbers)


# ### Exercise 2.8: 
# 
# Suppose arrays a and b are defined as follows:
# 
# from numpy import array
# 
# a = array([1,2,3,4],int)
# 
# b = array([2,4,6,8],int)
# 
# What will the computer print upon executing the following lines? (Try to work out the answer
# before trying it on the computer.)
# 

# In[59]:


from numpy import array
a = array([1,2,3,4],int)
b = array([2,4,6,8],int)

print(b/a+1)
print(b/(a+1))
print(1/a)


# ### Exercise 2.9: The Madelung constant
# 
# In condensed matter physics the Madelung constant gives the total electric potential felt by an atom in a solid. It depends on the charges on the other atoms nearby and their locations.
# 
# Consider for instance solid sodium chloride—table salt. The sodium chloride crystal has atoms arranged on a cubic lattice, but with alternating sodium and chlorine atoms, the sodium ones having a single positive charge +e and the chlorine ones a single negative charge −e, where e is
# the charge on the electron. If we label each position on the lattice by three integer coordinates (i, j, k), then the sodium atoms fall at positions where i + j + k is even, and the chlorine atoms at positions where i + j + k is odd.
# 
# Consider a sodium atom at the origin, i = j = k = 0, and let us calculate the Madelung constant. If the spacing of atoms on the lattice is a, then the distance from the origin to the atom at position (i, j, k) is 
# 
# $$\sqrt{((i\cdot a)^2) + ((j\cdot a)^2) + ((k\cdot a)^2)} = a\sqrt{i^2 + j^2 + k^2}$$
# 
# and the potential at the origin created by such an atom is
# 
# $$V(i,j,k) = \pm \frac{e}{4\pi\left(8.8541878128\times 10^{-12}\right)a\left(i^2 + j^2 + k^2\right)'}$$
# 
# with 8.8541878128×10^−12 being the permittivity of the vacuum and the sign of the expression depending on whether i + j + k is even or odd. The total potential felt by the sodium atom is then the sum of this quantity over all other atoms. Let us assume a cubic box around the sodium at the origin, with L atoms in all directions. 
# 
# where M is the Madelung constant, at least approximately—technically the Madelung constant is the value of M when L → ∞, but one can get a good approximation just by using a large value of L.
# Write a python program to calculate and print the Madelung constant for sodium chloride. Use as large a value of L as you can, while still having your program run in reasonable time say in a minute or less.

# In[62]:


import math

# Constants
e = 1.60217662e-19  # charge on an electron
epsilon_0 = 8.8541878128e-12  # permittivity of the vacuum
a = 5.64e-10  # lattice spacing
L = 8  # size of cubic box

# Function to calculate potential at the origin
def V(i, j, k):
    if (i + j + k) % 2 == 0:
        sign = 1
    else:
        sign = -1
    return sign * e / (4 * math.pi * epsilon_0 * a * math.sqrt(i**2 + j**2 + k**2))

# Madelung constant calculation
M = 0
for i in range(-L, L+1):
    for j in range(-L, L+1):
        for k in range(-L, L+1):
            if i == 0 and j == 0 and k == 0:
                continue
            M += V(i, j, k)

print("Madelung constant: ", M)


# Note that L=8 is used here as an example and it is good enough for reasonable time run as you mentioned in your question. The larger the value of L, the more accurate the Madelung constant will be, but it will also take longer for the program to run.

# ### Exercise 2.10: The semi-empirical mass formula
# 
# In nuclear physics, the semi-empirical mass formula is a formula for calculating the approximate nuclear binding energy B of an atomic nucleus with atomic number Z and mass number A is....
# 
# $$B = a_1A - a_2A^{2/3} - \frac{a_3(Z^2)}{A^{1/3}} - \frac{a_4(A-2Z)^2}{A} + \frac{a_5}{A^{1/2}}$$
# 
# 
# 
# $$
# a_{5} = \left\{ 
#     \begin{array} \ 
#           0   & \mbox{if A is odd , }\\ 
#          12.0 & \mbox{if A and Z are both even, }\\ 
#         −12.0 & \mbox{if A is even and Z is odd.} \\ 
#     \end{array} 
#     \right. 
# $$ 
# 
# a.) Write a program that takes as its input the values of A and Z, and prints out the binding
# energy for the corresponding atom. Use your program to find the binding energy of an
# atom with A = 58 and Z = 28. (Hint: The correct answer is around 490 MeV.)
# 
# b.) Modify your program to print out not the total binding energy B, but the binding energy per nucleon, which is B/A.
# 
# 
# c.) Now modify your program so that it takes as input just a single value of the atomic number Z and then goes through all values of A from A = Z to A = 3Z, to find the one that has the largest binding energy per nucleon. This is the most stable nucleus with the given atomic number. Have your program print out the value of A for this most stable nucleus and the value of the binding energy per nucleon.
# 
# d.) Modify your program again so that, instead of taking Z as input, it runs through all values of Z from 1 to 100 and prints out the most stable value of A for each one. At what value of Z does the maximum binding energy per nucleon occur? (The true answer, in real life, is Z = 28, which is nickel.)

# In[71]:


a1 = 15.8
a2 = 18.3
a3 = 0.714
a4 = 23.2

def binding_energy(A, Z):
    N = A - Z
    if A % 2 == 1:
        a5 = 0
    elif Z % 2 == 0:
        a5 = 12.0
    else:
        a5 = -12.0
    B = ((a1*A) - (a2*(A**(2/3))) - ((a3*(Z**2))/(A**(1/3))) - (a4*((A-(2*Z))**2)/A) + (a5/(A**(1/2))))
    return B

A = 58
Z = 28
binding_energy = binding_energy(A, Z)
print(binding_energy, "MeV")


# ### Exercise 2.10b
# 
# Modify your program to print out not the total binding energy B, but the binding energy per nucleon, which is B/A.

# In[72]:


a1 = 15.8
a2 = 18.3
a3 = 0.714
a4 = 23.2

def binding_energy(A, Z):
    N = A - Z
    if A % 2 == 1:
        a5 = 0
    elif Z % 2 == 0:
        a5 = 12.0
    else:
        a5 = -12.0
    B = ((a1*A) - (a2*(A**(2/3))) - ((a3*(Z**2))/(A**(1/3))) - (a4*((A-(2*Z))**2)/A) + (a5/(A**(1/2))))
    B_per_nucleon = B/A
    return B_per_nucleon

A = 58
Z = 28
binding_energy_per_nucleon = binding_energy(A, Z)
print(binding_energy_per_nucleon, "MeV/nucleon")


# ### Exercise 2.10c
# 
# Now modify your program so that it takes as input just a single value of the atomic number Z and then goes through all values of A from A = Z to A = 3Z, to find the one that has the largest binding energy per nucleon. This is the most stable nucleus with the given atomic number. Have your python program print out the value of A for this most stable nucleus and the value of the binding energy per nucleon.

# In[78]:


a1 = 15.8
a2 = 18.3
a3 = 0.714
a4 = 23.2

def binding_energy(A, Z):
    N = A - Z
    if A % 2 == 1:
        a5 = 0
    elif Z % 2 == 0:
        a5 = 12.0
    else:
        a5 = -12.0
    B = ((a1*A) - (a2*(A**(2/3))) - ((a3*(Z**2))/(A**(1/3))) - (a4*((A-(2*Z))**2)/A) + (a5/(A**(1/2))))
    B_per_nucleon = B/A
    return B_per_nucleon

Z = 28
max_B_per_nucleon = 0
max_A = 0
for A in range(Z, 3*Z):
    B_per_nucleon = binding_energy(A, Z)
    if B_per_nucleon > max_B_per_nucleon:
        max_B_per_nucleon = B_per_nucleon
        max_A = A

print("The most stable nucleus with atomic number", Z, "has mass number", max_A, "and binding energy per nucleon of", max_B_per_nucleon, "MeV/nucleon")


# In[80]:


a1 = 15.8
a2 = 18.3
a3 = 0.714
a4 = 23.2

def binding_energy(A, Z):
    N = A - Z
    if A % 2 == 1:
        a5 = 0
    elif Z % 2 == 0:
        a5 = 12.0
    else:
        a5 = -12.0
    B = ((a1*A) - (a2*(A**(2/3))) - ((a3*(Z**2))/(A**(1/3))) - (a4*((A-(2*Z))**2)/A) + (a5/(A**(1/2))))
    B_per_nucleon = B/A
    return B_per_nucleon

max_B_per_nucleon = 0
max_Z = 0
max_A = 0
for Z in range(1,101):
    for A in range(Z, 3*Z):
        B_per_nucleon = binding_energy(A, Z)
        if B_per_nucleon > max_B_per_nucleon:
            max_B_per_nucleon = B_per_nucleon
            max_Z = Z
            max_A = A

print("The most stable nucleus has atomic number", max_Z, "mass number", max_A, "and binding energy per nucleon of", max_B_per_nucleon, "MeV/nucleon")


# As you can see, the maximum binding energy per nucleon occurs at atomic number 28, which is Nickel. This is consistent with the expected result.
# 
# Please note that the SEMF does not work for all nuclei and its accuracy decrease for nuclei far from stability. It is only an approximation, and more complex models are needed for specific nuclei.

# ### Exercise 2.11: Binomial coefficients
# 
# Using this form for the binomial coefficient, write a user-defined function binomial(n,k) that calculates the binomial coefficient for given n and k. Make sure your function returns the answer in the form of an integer (not a float) and gives the correct value of 1 for the case where k = 0.

# In[82]:


def binomial(n, k):
    if k == 0:
        return 1
    else:
        return int(math.factorial(n) / (math.factorial(k) * math.factorial(n-k)))


# In[83]:


import math
def binomial(n, k):
    if k == 0:
        return 1
    else:
        return int(math.factorial(n) / (math.factorial(k) * math.factorial(n-k)))

for i in range(0, 21):
    for j in range(0, i+1):
        print(binomial(i, j), end = " ")
    print()


# In[84]:


import math

def binomial(n, k):
    return int(math.factorial(n) / (math.factorial(k) * math.factorial(n-k)))

# Total probability of getting heads 60 times
n = 100
k = 60
p = binomial(n, k) * (1/2)**k * (1/2)**(n-k)
print("Probability of getting heads 60 times:", p)

# Probability of getting heads 60 or more times
p_60_or_more = 0
for i in range(60, 101):
    p_60_or_more += binomial(n, i) * (1/2)**i * (1/2)**(n-i)
print("Probability of getting heads 60 or more times:", p_60_or_more)


# ### Exercise 2.12: Prime numbers
# 
# Develop a much faster program for prime numbers by making use of the following observations:
# 
# a) A number n is prime if it has no prime factors less than n. Hence we only need to check if it is divisible by other primes.
# 
# b) If a number n is non-prime, having a factor r, then n = rs, where s is also a factor. If r ≥ √n then n = rs ≥ √ns, which implies that s ≤ √n. In other words, any non-prime must have factors, and hence also prime factors, less than or equal to √n. Thus to determine if a number is prime we have to check its prime factors only up to and including √n—if there  are none then the number is prime.
# 
# c) If we find even a single prime factor less than √n then we know that the number is non-prime, and hence there is no need to check any further—we can abandon this number and move on to something else.
# 
# 
# Write a Python program that finds all the primes up to ten thousand. Create a list to store the primes, which starts out with just the one prime number 2 in it. Then for each number n from 3 to 10 000 check whether the number is divisible by any of the primes in the list up to and including √n. As soon as you find a single prime factor you can stop checking the rest of them—you know n is not a prime. If you find no prime factors √n or less then n is prime and you should add it to the list. You can print out the list all in one go at the end of the program,
# or you can print out the individual numbers as you find them.

# In[93]:


import math

# Initialize the list of primes with the first prime, 2
primes = [2]

# Iterate through each number from 3 to 10,000
for n in range(3, 10000):
    # Assume the number is prime
    is_prime = True

    # Iterate through each prime in the list up to and including sqrt(n)
    for prime in primes:
        if prime > math.sqrt(n):
            break
        # If the number is divisible by a prime in the list, it is not prime
        if n % prime == 0:
            is_prime = False
            break
    # If the number is prime, add it to the list
    if is_prime:
        primes.append(n)

# Print out the list of primes
print(primes)


# ### Exercise 2.13: Recursion

# In[85]:


def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)


# In[89]:


catalan_numbers = []

def catalan(n):
    if n == 0:
        return 1
    else:
        return int((4*n-2)/(n+1)*catalan(n-1))

n = 0
while catalan(n) <= 1000000000:
    catalan_numbers.append(catalan(n))
    n += 1

print(catalan_numbers)
print(catalan(100))


# In[91]:


def gcd(m, n):
    if n == 0:
        return m
    else:
        return gcd(n, m%n)

print(gcd(108, 192))


# ### Practising with Latex

# $\$ $ A_{n} = 10\sin\left( n\theta \right) $ \$ $ 
# 
# $ \left( \begin{array} &n \\ 0 \end{array} \right) = 1$

# $$
# A_{j} = \left\{ 
#     \begin{array} \ 
#           0   & \mbox{if j }\lt\mbox{ 10,} \\ 
#          12.0 & \mbox{ if j }\geq \mbox{ 10 and }\leq \mbox{ 20,} \\ 
#         −12.0 & \mbox{if j }\geq \mbox{ 20.} \\ 
#     \end{array} 
#     \right. 
# $$ 

# \begin{align} \ 
# Y\left( x, t\right) = A*cos(kx-\omega t) \ 
# \end{align} 
#  
# \begin{align} \ 
# A_{j} = \left\{ 
#     \begin{array} \ 
#           0   & \mbox{if j }\lt\mbox{ 10,} \\ 
#          12.0 & \mbox{ if j }\geq \mbox{ 10 and }\leq \mbox{ 20,} \\ 
#         −12.0 & \mbox{if j }\geq \mbox{ 20.} \\ 
#     \end{array} 
#     \right. 
# \end{align} 

# $$
# sign(x) = \left\{
#     \begin{array}\\
#         1 & \mbox{if } \ x \in \mathbf{N}^* \\
#         0 & \mbox{if } \ x = 0 \\
#         -1 & \mbox{else.}
#     \end{array}
# \right.
# $$
# 
# \\
# 
# $$
#  \left.
#     \begin{array} \\
#         \alpha^2 = \sqrt5 \\
#         \alpha \geq 0 
#     \end{array}
# \right \}=\alpha = 5 
# $$

# To insert a mathematical formula we use the dollar symbol $, as follows:
# 
# Euler's identity: $ e^{i \pi} + 1 = 0 $
# 
# To isolate and center the formulas and enter in math display mode, we use 2 dollars symbol:
# 
# $$
# ...
# $$
# 
# 
# 
# Euler's identity: $$ e^{i \pi} + 1 = 0 $$
# 
# 

# To add little spacing in math mode use \,
# 
# To add a new line when in math mode use \\
# 
# To display fraction use \frac{arg 1}{arg 2}
# 
# For power (superscripts text) use ^{}
# 
# For indices (subscripts) use _{}
# 
# For roots use \sqrt[n] {arg}
# 
# The [n]is optional.
# 

# Given : $\pi = 3.14$ , $\alpha = \frac{3\pi}{4}\, rad$
# $$
# \omega = 2\pi f \\
# f = \frac{c}{\lambda}\\
# \lambda_0=\theta^2+\delta\\
# \Delta\lambda = \frac{1}{\lambda^2}
# $$

# ## Chapter 3 Examples

# In[7]:


from pylab import plot,show
y = [ 1.0, 2.4, 1.7, 0.3, 0.6, 1.8]
plot (y)
show()


# In[9]:


from pylab import plot,show
x = [ 0.5, 1.0, 2.0, 4.0, 7.0, 10.0 ] 
y = [ 1.0, 2.4, 1.7, 0.3, 0.6, 1.8] 
plot(x,y)
show()


# In[12]:


from pylab import plot, show 
from numpy import linspace, sin
x =linspace(0,10,100) 
y = sin(x)
plot(x,y)
show()


# In[17]:


from numpy import loadtxt 
from pylab import plot, show

#data = loadtxt("values.txt",float)
#x = data[:,0] 
#y = data[:,1]
#plot(x,y) 
#show ()


# In[ ]:




