import random 
import math
from operator import add 
from functools import reduce
import matplotlib.pyplot as plt 
#Defining constants 
a=5.3*(10**-10)
sigma=3.4*(10**-10)
epsilon=10*(10**3)
alpha=0.2
# Kb=1.38*(10**(−10))
Kb=1.3*(10**-23)
T=5
e=2.71828
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
energy_original=[]
average_original=[]
energy_prime=[]
acceptance=[]
x=[]
y=[]
#MC Method
for ri in particle_position:
    for j in particle_position:
        if ri!=j:
            rij=calculate_distance(ri,j)
            if rij<(3*sigma):
                energy_original.append(calc_energy(rij))
    Ui=average_energy(energy_original)
    average_original.append(Ui)
    for imc in (0,1):
        ri_prime=r_prime_calc(ri)
        # for j in particle_position and j!=ri:
        for j in particle_position:
                if j!=ri:
                    rij_prime=calculate_distance(ri_prime,j)
                    if rij_prime<(3*sigma):
                        energy_prime.append(calc_energy(rij_prime))
        Ui_prime=average_energy(energy_prime)
        average_prime.append(Ui_prime)
        if Ui_prime<Ui:
            acceptance.append(1)
        else:
            delta=Ui_prime-Ui
            P = e**(-delta/(Kb*T))
            G = random.random()
            if G<P:
                acceptance.append(1)
            else:
                acceptance.append(0)
    print(acceptance)
    break

            



    # for i in range(0,1):
    #     r_prime=r_prime_calc(r)
    #     for j in particle_position:
    #         if j!=r:
    #             rij=calculate_distance(r,j)
    #             rij_prime=calculate_distance(r_prime,j)
    #             # print(rij_prime,rij)
    #             if float(rij)<(3*sigma):
    #                 phi_rij=calc_energy(rij)
    #                 x.append(rij)
    #                 y.append(phi_rij)
    #                 energy_original.append(phi_rij)
    #             # if float(rij_prime)<(3*sigma): #Bug -> must be there actually 
    #             phi_rij_prime=calc_energy(rij_prime)
    #             energy_prime.append(phi_rij_prime)
            

    #             if phi_rij>phi_rij_prime:
    #                 acceptance.append(1)
    #             else:
    #                 delta=phi_rij_prime-phi_rij
    #                 P = e**(-delta/(Kb*T))
    #                 G = random.random()
    #                 if G<P:
    #                     acceptance.append(1)
    #                 else:
    #                     acceptance.append(0) 
                    
    #     avg_original=average_energy(energy_original)
    #     avg_prime=average_energy(energy_prime)
    #     change.append([r,r_prime,avg_original,avg_prime])
# print(len(average_original))
# probability=0
# for i in acceptance:
#     if i==1:
#         probability=probability+1
# probability=probability/len(acceptance)
# print("The probability of acceptance is:",probability,"%")

