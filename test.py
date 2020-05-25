import random 
import math
from operator import add 
from functools import reduce
#Defining constants 
a=5.3*(10**-10)
sigma=3.4*(10**-10)
epsilon=10*(10**3)
alpha=0.2
#Generating N=500 data using fcc base
def create_N_particles( n ):
    x=(n/4)**(1/3)
    particle_position=[]
    for i in range(0,5):
        for j in range(0,5):
            for k in range (0,5):
                particle_position.append([0+i,0+j,0+k])
                particle_position.append([0.5+i,0.5+j,0+k])
                particle_position.append([0.5+i,0+j,0.5+k])
                particle_position.append([0+i,0.5+j,0.5+k])
    for i in particle_position:
        for j in range(0,3):
            i[j]=a*i[j]
    return particle_position

def calculate_distance(ri,rj):
    a=(ri[0]-rj[0])**2
    b=(ri[1]-rj[1])**2
    c=(ri[2]-rj[2])**2
    rij=math.sqrt(a+b+c)
    return rij

def calc_energy(r):
    phi_rij=4*epsilon*(((sigma/r)**12)-((sigma/r)**6))
    return phi_rij

def average_energy(lst): 
    return reduce(lambda a, b: a + b, lst) / len(lst) 

def r_prime_calc(r):
    n1=((random.random()*2)-1)*alpha
    n2=((random.random()*2)-1)*alpha
    n3=((random.random()*2)-1)*alpha
    random_nums=[n1,n2,n3]
    r_prime=list(map(add, r, random_nums))
    return r_prime


particle_position=create_N_particles(500)
change=[]
energy=[]

# # for first particle 
# r=particle_position[0]
# for i in range(0,5):
#     r_prime=r_prime_calc(r)
#     for j in particle_position:
#         if j!=r:
#             rij=calculate_distance(r,j)
#             if float(rij)<(3*sigma):
#                 phi_rij=calc_energy(rij)
#                 energy.append(phi_rij)
#     i=i+1

#     avg=average_energy(energy)
#     change.append([r,r_prime,avg])
# print(change)
                

for r in particle_position:
    for i in range(0,5):
        r_prime=r_prime_calc(r)
        for j in particle_position:
            if j!=r:
                rij=calculate_distance(r,j)
                if float(rij)<(3*sigma):
                    phi_rij=calc_energy(rij)
                    energy.append(phi_rij)
        avg=average_energy(energy)
        change.append([r,r_prime,avg])
print(len(change))


