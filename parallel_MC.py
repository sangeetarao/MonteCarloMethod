import random 
import math
from operator import add 
import matplotlib.pyplot as plt 
import time
start=time.time()

#Defining constants and initializing 
a=5.3*(10**-10)
L=5*a
sigma=3.4*(10**-10)
epsilon=10**3
epsilon_4=4*epsilon
alpha=0.005*(10**-10) #should be around 0.3 A
R=8.3145
T=100

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
#Defining variables
particle_position=create_N_particles(500)
energy_original=[]
acceptance=0
Ui=[]
mc_cycle=60000
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
print("Using neighbours methods")
print("Energy of base fcc configuration is: ", (U/500)/1000,"KJ/mol")
Uri=U

for imc in range(0,mc_cycle):
    print("MC Cycle:",imc)
    print("\tEnergy of OG configuration is: ", (Uri/500)/1000,"KJ/mol")
    pp_prime=[]
    Uri_list=[]
    Uri_prime_list=[]
    if imc%10==0:
        neighbours=[]
        for ri in particle_position:
            n_ri=[]
            for j in particle_position:
                if ri!=j:
                    rij=calculate_distance(ri,j)
                    if rij<((3*sigma)+(1*(10**-10))):
                        n_ri.append(particle_position.index(j))
            neighbours.append(n_ri)
            # print(neighbours[0])

    for ri in particle_position:
        pp_prime.append(r_prime_calc(ri))        
    for ri in pp_prime:
        energy_prime=[]
        n_index=pp_prime.index(ri)
        for j_index in neighbours[n_index]:
                j=pp_prime[j_index]
                if ri!=j:
                    rij=calculate_distance(ri,j)
                    if rij<(3*sigma):
                        energy_prime.append(calc_energy(rij))
        Uri_prime_list.append(sum(energy_prime))
    # Uri=sum(Uri_list)
    Uri_p=sum(Uri_prime_list)
    delta=Uri_p-Uri
    if Uri_p<Uri:
        particle_position=pp_prime.copy()
        acceptance=acceptance+1
        Uri=Uri_p
        pass 
    else:
        P=math.exp((-delta)/(R*T))
        G=random.random()
        if G<P:
            particle_position=pp_prime.copy()
            acceptance=acceptance+1
            Uri=Uri_p
        else:
            pass
    avg_accept=acceptance/(imc+1)
    print("\tEnergy of PR configuration is: ", (Uri_p/500)/1000,"KJ/mol")
    print("\tAvg acceptance:",avg_accept*100,"i.e",acceptance,"/",imc+1)
    x.append(imc)
    y.append((Uri/500))
end=time.time()
print("Total time taken:",end-start)
print("Particle position")
print(particle_position)
print("Energy list: (for final graph)")
print(y)
plt.plot(x,y)
plt.title("30-60k with 0.005 alpha")
plt.xlabel("MC Cylces")
plt.ylabel("Average energy per particle")
plt.savefig("30-60k.png")




