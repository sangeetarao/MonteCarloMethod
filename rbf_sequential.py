import random 
import math
from operator import add 
import matplotlib.pyplot as plt 
import time
start=time.time()

#Generating N=500 data using fcc lattice
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
	dx=(to_make_bulk(ri[0],rj[0]))**2
	dy=(to_make_bulk(ri[1],rj[1]))**2
	dz=(to_make_bulk(ri[2],rj[2]))**2
	rij=math.sqrt(dx+dy+dz)
	return rij

def to_make_bulk(i,j):
	dx=i-j 
	if dx>L/2:
		dx=dx-L 
	if dx<(-L/2):
		dx=dx+L 
	return dx

def calc_energy(r):
	constant=(sigma/r)**6
	phi_rij=epsilon_4*((constant**2)-constant)
	return phi_rij

def r_prime_calc(r):
    n1=((random.random()*2)-1)*alpha
    n2=((random.random()*2)-1)*alpha
    n3=((random.random()*2)-1)*alpha
    random_nums=[n1,n2,n3]
    r_prime=pbc(list(map(add, r, random_nums)))
    return r_prime

def pbc(r):
	if r[0]>L:
		r[0]=r[0]-L 
	elif r[0]<0:
		r[0]=r[0]+L
	if r[1]>L:
		r[1]=r[1]-L 
	elif r[1]<0:
		r[1]=r[1]+L
	if r[2]>L:
		r[2]=r[2]-L 
	elif r[2]<0:
		r[2]=r[2]+L
	return r


#Defining constants and initializing 
a=5.3*(10**-10)
L=5*a
sigma=3.4*(10**-10)
epsilon=10**3
epsilon_4=4*epsilon
alpha=0.3*(10**-10) #should be around 0.3 A
R=8.3145
T=100
particle_position=create_N_particles(500)
energy_original=[]
acceptance=0
Ui=[]
mc_cycle=5
x=[]
y=[]
#MC METHOD

#Energy of system for base configuration
for ri in particle_position:
	energy_original=[]
	for j in particle_position:
			if ri!=j:
				rij=calculate_distance(ri,j)
				if rij<(3*sigma):
					energy_original.append(calc_energy(rij))
	Ui.append(sum(energy_original))
U=sum(Ui)
U_run=U
print("Energy of base fcc configuration is: ", (U/500)/1000,"KJ/mol")

#Code with neighbours 


#Starting MC Cycles
for imc in range(0,mc_cycle):
	#To calculate neighbours for every 10 cycles
	if imc%10==0:
		neighbours=[]
		for ri in particle_position:
			n_for_ri=[]
			for j in particle_position:
				if ri!=j:
					rij=calculate_distance(ri,j)
					if rij<(3*sigma)+(1*(10**-10)):
						n_for_ri.append(particle_position.index(j))
			neighbours.append(n_for_ri)
		print("--------------------------------------------------------------------------------------------------------------------------------")

	#To calculate energy of system before moving particle
	energy_of_system=[]
	for ri in particle_position:
		eng_og=[]
		eng_prime=[]
		n_index=particle_position.index(ri)
		for j_index in neighbours[n_index]:
			j=particle_position[j_index]
			if ri!=j:
				rij=calculate_distance(ri,j)
				if rij<(3*sigma):
					eng_og.append(calc_energy(rij))
		Uri=(sum(eng_og))
		energy_of_system.append(Uri)
		#Generating r prime
		ri_prime=r_prime_calc(ri)
		for j_index in neighbours[n_index]:
			j=particle_position[j_index]
			if ri!=j:
				rij_prime=calculate_distance(ri_prime,j)
				if rij_prime<(3*sigma):
					eng_prime.append(calc_energy(rij_prime))
		#energy of system with r prime (after move)
		Uri_prime=(sum(eng_prime))
		delta=Uri_prime-Uri
		#if move is accepted then change Urun by adding delta
		if Uri_prime<Uri:
			acceptance=acceptance+1
			index=particle_position.index(ri)
			particle_position[index]=ri_prime
			U_run+=(delta)
		else:
			P = math.exp(-delta/(R*T))
			G = (random.random())
			if G<P:
				acceptance=acceptance+1
				index=particle_position.index(ri)
				particle_position[index]=ri_prime
				U_run+=(delta)
			#Reject if move is not accepted
			else:
				pass.
	#to calculate new alpha based on average acceptance rate 
	x.append(imc)
	y.append(U_run/500)
	alpha_test=acceptance/((imc+1)*500)
	#to reset Urun for every 20 cycles
	# if (imc!=0 and imc%20==0):
	# 	U_run=sum(energy_of_system)
	if (imc!=0 and imc%10==0):
		# print("Energy of system(wrt delta)",(U_run/500)/1000,"Energy of system by calculation",(sum(energy_of_system)/500)/1000)
		if alpha_test<0.5:
			alpha=alpha-(0.01*(10**-10))
			print("Changed alpha to:",alpha,"Energy of system(wrt delta)",(U_run/500)/1000,"Energy of system by calculation",(sum(energy_of_system)/500)/1000)
		elif alpha_test>0.5:
			alpha=alpha+(0.01*(10**-10))
			print("Changed alpha:",alpha,"Energy of system(wrt delta)",(U_run/500)/1000,"Energy of system by calculation",(sum(energy_of_system)/500)/1000)
	print("MC Cycle:",imc)
	print("\tTotal Energy (Average for ith particle):",(U_run/1000)/500,"KJ/mol")
	print(len(neighbours[0]))
end=time.time()
print("\nTotal acceptance is:",alpha_test,"Final alpha value:",alpha,"Time taken: ",end-start)
#To plot energy vs mc cycle 
plt.plot(x,y)
plt.show()


