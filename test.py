import random 
import math
#Defining constants 
a=5.3*(10**-10)
sigma=3.4*(10**-10)
epsilon=10*(10**3)
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
    rij=((ri[0]-rj[0])**2+(ri[1]-rj[1])**2+(ri[2]-rj[2]**2))**0.5
    return rij

def calc_energy(r):
    phi_rij=4*epsilon*(((sigma/r)**12)-((sigma/r)**6))
    return phi_rij

def average_energy(lst): 
    return reduce(lambda a, b: a + b, lst) / len(lst) 



particle_position=create_N_particles(500)
print(len(particle_position))
sample_list=[[0,0,0],[0,1,0],[0,0,1],[1,0,0]]
changes=[]
energy=[]
alpha=0.2
# Starting MC Method trial 
ri=sample_list[0]
rj=sample_list[1]
rij=calculate_distance(ri,rj)
energy.append(calc_energy(rij))

print(energy)
n1=(random.random()*2)-1
n2=(random.random()*2)-1
n3=(random.random()*2)-1

print(n1,n2,n3)
