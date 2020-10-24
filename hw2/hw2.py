import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt 
import math
import random

# section1

#defining variables, all units in CGS, 
#100 pc = 3.086e+20 cm
D = 3.086e+20 #cm
n = 1 #cm-3
# using column density (N) = n*D
N = D*n
print("column density = ",N, "cm^-2")
#using Optical depth (tau) = N*cross section (sigma)
tautot = [10**-3,1,10**3]
Sig = [a/N for a in tautot]
print ("Sigma =", Sig, "for: Tau = 10^-3,","Tau =1,", 'Tau=10^3')

#section 2

def steps(D,ds):
    steparray = np.linspace(0,D,ds)
    print (steparray)

 #defining the array holding the Specific Intensity
array2 =[]
def Iv(Sv,I0,Tau):
    IvTau = [I0*math.e**(-a)+(Sv*math.e**(-a))*(math.e**(a)-1) for a in Tau]
    return IvTau
    print(IvTau)

#testing function
t=np.linspace(0,100,1000)
IvTest = Iv(1,5,t)
for i in range(1000):
    array2.append(Iv(1,5,t))
array2

#section 3
#defining Sigma\nu as a gaussian function (y=sigma0*e^(-x-freq0)/(2*stdev)^2)) for sigma0 occuring at some frequency freq0

F = np.linspace(0,100,100)
def sigmanu(sigma0,mu,f,f0):
    y=[sigma0*math.e**(-(a-f0)**2/((2*mu)**2)) for a in f]
    return y
    print(y)

#testing function
sigma1 = sigmanu(10e-3,2,F,50)

#plot of sigma0 = 3.240440699935191e-24
plt.plot(F, sigmanu(3.240440699935191e-24,2,F,50))
plt.title('sigma0 = 3.240440699935191e-24')
plt.xlabel('Frequency')
plt.ylabel('Sigma')
plt.show()

#plot of sigma0 = 3.240440699935191e-21
plt.plot(F, sigmanu(3.240440699935191e-21,2,F,50))
plt.title('sigma0 = 3.240440699935191e-21')
plt.xlabel('Frequency')
plt.ylabel('Sigma')
plt.show()

#plot of sigma0 = 3.2404406999351913e-18
plt.plot(F, sigmanu(3.2404406999351913e-18,2, F, 50))
plt.title('sigma0 = 3.2404406999351913e-18')
plt.xlabel('Frequency')
plt.ylabel('Sigma')
plt.show()

#section 4

#a
print("Iv(0)=0,tau(D)<1")
Taua = [N*a for a in sigmanu(3.240440699935191e-24,2,F,50)]
I0a=0
Sva=1
#Iv(Sva,I0a,Taua)
plt.plot(F,Iv(Sva,I0a,Taua))
plt.title("Iv(0)=0,tau(D)<1")
plt.xlabel('Frequency')
plt.ylabel('Iv')
plt.show()

#b
print("Iv(0)>Sv,tau(D)<1")
Taub = [N*a for a in sigmanu(3.240440699935191e-24,2,F,50)]
I0b=1
Svb=0.5
Iv(Svb,I0b,Taub)
plt.plot(F,Iv(Svb,I0b,Taub))
plt.title("Iv(0)>Sv,tau(D)<1")
plt.xlabel('Frequency')
plt.ylabel('Iv')
plt.show()

#c
print("Iv(0)<Sv,tau(D)<1")
Tauc = [N*a for a in sigmanu(3.240440699935191e-24,2,F,50)]
I0c=0.5
Svc=1
Iv(Svc,I0c,Tauc)
plt.plot(F,Iv(Svc,I0c,Tauc))
plt.title("Iv(0)<Sv,tau(D)<1")
plt.xlabel('Frequency')
plt.ylabel('Iv')
plt.show()

#d
print("tau(D)>>1")
Taud = [100 for a in sigmanu(3.2404406999351913e-18,2,F,50)]
I0d=0.5
Svd=1
Iv(Svd,I0d,Taud)
plt.plot(F,Iv(Svd,I0d,Taud))
plt.title("tau(D)>>1")
plt.xlabel('Frequency')
plt.ylabel('Iv')
plt.show()

#e
print("Iv(0)<Sv,tau(D)<1 while tau0(D)>1")
Taue = [N*a for a in sigmanu(6.240440699935191e-21,2,F,50)]
I0e=0.5
Sve=1
Iv(Sve,I0e,Taue)
plt.plot(F,Iv(Sve,I0e,Taue))
plt.title("Iv(0)<Sv,tau(D)<1 while tau0(D)>1")
plt.xlabel('Frequency')
plt.ylabel('Iv')
plt.show()

#f
print("Iv(0)>Sv,tau(D)<1 while tau0(D)<1")
Tauf = [N*a for a in sigmanu(6.240440699935191e-21,2,F,50)]
I0f=1
Svf=0.5
Iv(Svf,I0f,Tauf)
#plt.plot(F,Tauf)
plt.plot(F,Iv(Svf,I0f,Tauf))
plt.title("Iv(0)>Sv,tau(D)<1 while tau0(D)>1")
plt.xlabel('Frequency')
plt.ylabel('Iv')
plt.show()

