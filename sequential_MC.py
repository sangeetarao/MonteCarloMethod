import random 
import math
from operator import add 
from functools import reduce
import matplotlib.pyplot as plt 

#Generating N=500 data using fcc base
def create_N_particles( n ):
    x=int((n/4)**(1/3))
    particle_position=[]
    for i in range(0,x):
        for j in range(0,x):
            for k in range (0,x):
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

def r_prime_calc(r):
    n1=((random.random()*2)-1)*alpha
    n2=((random.random()*2)-1)*alpha
    n3=((random.random()*2)-1)*alpha
    random_nums=[n1,n2,n3]
    r_prime=list(map(add, r, random_nums))
    return r_prime
def calc_acceptance(lst):
	acceptance_rate=sum(lst)/len(lst)
	return acceptance_rate




#Defining constants and initializing 
a=5.3*(10**-10)
sigma=3.4*(10**-10)
epsilon=10**3
alpha=0.3*(10**-10) #should be around 0.3 A
R=8.3145
T=100
e=2.718281828
particle_position=create_N_particles(500)
energy_original=[]
acceptance=[]
Ui=[]
configurations=[]
mc_cycle=100
#MC METHOD
test=particle_position.copy()
#Energy of system for base configuration
for ri in test:
    energy_original=[]
    for j in test:
        if ri!=j:
            rij=calculate_distance(ri,j)
            if rij<(3*sigma):
                energy_original.append(calc_energy(rij))
    Ui.append(sum(energy_original))
U=sum(Ui)/500
U_run=U
print("Energy of base configuration is: ", U ,"J/mol")

eng_og=[]
eng_prime=[]
check=0
#Starting MC Cycles
for imc in range(0,mc_cycle):
	for ri in test:
		eng_og=[]
		eng_prime=[]
		for j in test:
			if ri!=j:
				rij=calculate_distance(ri,j)
				if rij<(3*sigma):
				    eng_og.append(calc_energy(rij))
		Uri=sum(eng_og)
		ri_prime=r_prime_calc(ri)
		for j in test:
			if j!=ri:
				rij_prime=calculate_distance(ri_prime,j)
				if rij_prime<(3*sigma):
				    eng_prime.append(calc_energy(rij_prime))
		Uri_prime=sum(eng_prime)
		delta=Uri_prime-Uri
		# print(delta)
		if Uri_prime<Uri:
			acceptance.append(1)
			check=1
		else:
			P = math.exp(-delta/(R*T))
			G = (random.random())*2-1
			if G<P:
				acceptance.append(1)
				check=1
			else:
				acceptance.append(0)
				check=0
		if check==1:
			index=test.index(ri)
			test[index]=ri_prime
			U_run=U_run+delta
	if (imc%10)==0:
		alpha_test=calc_acceptance(acceptance)
		if alpha_test<0.5:
			alpha=alpha-(0.01*(10**-10))
			print("Changed alpha to: ",alpha)
		elif alpha_test>0.5:
			alpha=alpha+(0.01*(10**-10))
			print("Changed alpha to: ",alpha)
	configurations.append(test)
	print("MC Cycle: ",imc+1," Energy:", U_run," Acceptance:",calc_acceptance(acceptance))


print("Total acceptance is:",calc_acceptance(acceptance),"Final alpha value:",alpha)

        



