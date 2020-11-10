import numpy as np
from numpy import *
import scipy as scipy
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt 
import math
import random

# all units are in CGS
a0 = 5.2918*(10**-9)
e =  4.80320425*(10**-10)
Z = 2
me = 9.109*(10**-24)
IonCharge = Z*e
xi_ion = 0
yi_ion = 0
xi_e = a0*500
yi_e = a0*500
vi_e = 10**7 
t_i = 0

#Setting axes and timescale
x_ion = np.linspace(500,0.1,5000)
y_ion = np.linspace(500,0.1,5000)
x_e = np.linspace(xi_e,0.1,5000)
y_e = np.linspace(yi_e,0.1,5000)
r_e = [math.sqrt(a**2 +b**2) for a,b in zip(x_e,y_e)]
time = np.linspace(1e-22,1e-10,5000)


def calcs(x_i, v_i, ti, deltat, e, IonCharge, me):
        x_i = x_i + x
        t_i = t_i + deltat
        f = e*IonCharge/x_i**2
        a_x = f/me
        v = a_x*t_i
        x = v*t_i - .5*a_x*t_i**2
        return a_x
        return v
        return x
        print(a_x)

#using Cuolomb's law F= m*a = abs(q1*q2)/r^2 in CGS
F_ix = e*IonCharge/(xi_e**2)
a_ix= F_ix/me
v_2 = vi_e + a_ix*.005
x_2 = xi_e + v_2*.005 -.5*a_ix*.005**2
f_2 = e*IonCharge/x_2
print(f_2)

def position_x(v_i, t, e, IonCharge, me, x_i):
    N = t.size - 1
    dt = t[1] - t[0]  
    u = np.zeros(N+1)
    v = np.zeros(N+1)
    f = np.zeros(N+1)
    a = np.zeros(N+1)
    u[0] = x_i
    f[0] = e*IonCharge/(u[0]**2)
    a[0] = f[0]/me
    v[0] = v_i
    v[1] = v[0] + a[0]*t[1]
    u[1] = u[0] + v[1]*t[1] - .5*a[1]*t[1]**2
    f[1] = e*IonCharge/(u[1]**2)
    a[1] = f[1]/me
    
    for n in range(1,N):
        f[n+1] = e*IonCharge/(u[n]**2)
        a[n+1] = f[n+1]/me
        v[n+1] = v[n] +a[n+1]*t[n+1]
        u[n+1] = u[n] + v[n+1]*t[n+1] - .5*a[n+1]*t[n+1]**2
        
    return u
    #return f
    #return a
    #return v
    print (u)
    #print(f)
    #print(a)
    #print(v)


def force_x(v_i, t, e, IonCharge, me, x_i):
    N = t.size - 1
    dt = t[1] - t[0]  
    u = np.zeros(N+1)
    v = np.zeros(N+1)
    f = np.zeros(N+1)
    a = np.zeros(N+1)
    u[0] = x_i
    f[0] = e*IonCharge/(u[0]**2)
    a[0] = f[0]/me
    v[0] = v_i
    v[1] = v[0] + a[1]*t[1]
    u[1] = u[0] + v[1]*t[1] - .5*a[1]*t[1]**2
    f[1] = e*IonCharge/(u[1]**2)
    a[1] = f[1]/me
    
    for n in range(1,N):
        f[n+1] = e*IonCharge/(u[n]**2)
        a[n+1] = f[n+1]/me
        v[n+1] = v[n] +a[n+1]*t[n+1]
        u[n+1] = u[n] + v[n+1]*t[n+1] - .5*a[n+1]*t[n+1]**2
        
    #return u
    return f
    #return a
    #return v
    #print (u)
    print(f)
    #print(a)
    #print(v)
def acc_x(v_i, t, e, IonCharge, me, x_i):
    N = t.size - 1
    dt = t[1] - t[0]  
    u = np.zeros(N+1)
    v = np.zeros(N+1)
    f = np.zeros(N+1)
    a = np.zeros(N+1)
    u[0] = x_i
    f[0] = e*IonCharge/(u[0]**2)
    a[0] = f[0]/me
    v[0] = v_i
    v[1] = v[0] + a[1]*t[1]
    u[1] = u[0] + v[1]*t[1] - .5*a[1]*t[1]**2
    f[1] = e*IonCharge/(u[1]**2)
    a[1] = f[1]/me
    
    for n in range(1,N):
        f[n+1] = e*IonCharge/(u[n]**2)
        a[n+1] = f[n+1]/me
        v[n+1] = v[n] +a[n+1]*t[n+1]
        u[n+1] = u[n] + v[n+1]*t[n+1] - .5*a[n+1]*t[n+1]**2
        
    #return u
    #return f
    return a
    #return v
    #print (u)
    #print(f)
    print(a)
    #print(v)
def vel_x(v_i, t, e, IonCharge, me, x_i):
    N = t.size - 1
    dt = t[1] - t[0]  
    u = np.zeros(N+1)
    v = np.zeros(N+1)
    f = np.zeros(N+1)
    a = np.zeros(N+1)
    u[0] = x_i
    f[0] = e*IonCharge/(u[0]**2)
    a[0] = f[0]/me
    v[0] = v_i
    v[1] = v[0] + a[1]*t[1]
    u[1] = u[0] + v[1]*t[1] - .5*a[1]*t[1]**2
    f[1] = e*IonCharge/(u[1]**2)
    a[1] = f[1]/me
    
    for n in range(1,N):
        f[n+1] = e*IonCharge/(u[n]**2)
        a[n+1] = f[n+1]/me
        v[n+1] = v[n] +a[n+1]*t[n+1]
        u[n+1] = u[n] + v[n+1]*t[n+1] - .5*a[n+1]*t[n+1]**2
        
    #return u
    #return f
    #return a
    return v
    #print (u)
    #print(f)
    #print(a)
    print(v)

x = position_x(vi_e, time, e, IonCharge, me, xi_e)
print(x)
f = force_x(vi_e, time, e, IonCharge, me, xi_e)
print(f)
a = acc_x(vi_e, time, e, IonCharge, me, xi_e)
print(a)
v = vel_x(vi_e, time, e, IonCharge, me, xi_e)
print(v)

plt.plot(time,x)
plt.xlabel("time")
plt.ylabel("x position")
plt.title("X Postion verus time")
plt.show()

plt.plot(time,a)
plt.xlabel("time")
plt.ylabel("x acceleration")
plt.title("X acceleration verus time")
plt.show()

plt.plot(time,v)
plt.xlabel("time")
plt.ylabel("x velocity")
plt.title("X velocity verus time")
plt.show()

def position_y(v_i, t, e, IonCharge, me, y_i):
    N = t.size - 1
    dt = t[1] - t[0]  
    u = np.zeros(N+1)
    v = np.zeros(N+1)
    f = np.zeros(N+1)
    a = np.zeros(N+1)
    u[0] = y_i
    f[0] = e*IonCharge/(u[0]**2)
    a[0] = f[0]/me
    v[0] = v_i
    v[1] = v[0] + a[0]*t[1]
    u[1] = u[0] + v[1]*t[1] - .5*a[1]*t[1]**2
    f[1] = e*IonCharge/(u[1]**2)
    a[1] = f[1]/me
    
    for n in range(1,N):
        f[n+1] = e*IonCharge/(u[n]**2)
        a[n+1] = f[n+1]/me
        v[n+1] = v[n] +a[n+1]*t[n+1]
        u[n+1] = u[n] + v[n+1]*t[n+1] - .5*a[n+1]*t[n+1]**2
        
    return u
    #return f
    #return a
    #return v
    print (u)
    #print(f)
    #print(a)
    #print(v)
def force_y(v_i, t, e, IonCharge, me, y_i):
    N = t.size - 1
    dt = t[1] - t[0]  
    u = np.zeros(N+1)
    v = np.zeros(N+1)
    f = np.zeros(N+1)
    a = np.zeros(N+1)
    u[0] = y_i
    f[0] = e*IonCharge/(u[0]**2)
    a[0] = f[0]/me
    v[0] = v_i
    v[1] = v[0] + a[1]*t[1]
    u[1] = u[0] + v[1]*t[1] - .5*a[1]*t[1]**2
    f[1] = e*IonCharge/(u[1]**2)
    a[1] = f[1]/me
    
    for n in range(1,N):
        f[n+1] = e*IonCharge/(u[n]**2)
        a[n+1] = f[n+1]/me
        v[n+1] = v[n] +a[n+1]*t[n+1]
        u[n+1] = u[n] + v[n+1]*t[n+1] - .5*a[n+1]*t[n+1]**2
        
    #return u
    return f
    #return a
    #return v
    #print (u)
    print(f)
    #print(a)
    #print(v)
def acc_y(v_i, t, e, IonCharge, me, y_i):
    N = t.size - 1
    dt = t[1] - t[0]  
    u = np.zeros(N+1)
    v = np.zeros(N+1)
    f = np.zeros(N+1)
    a = np.zeros(N+1)
    u[0] = y_i
    f[0] = e*IonCharge/(u[0]**2)
    a[0] = f[0]/me
    v[0] = v_i
    v[1] = v[0] + a[1]*t[1]
    u[1] = u[0] + v[1]*t[1] - .5*a[1]*t[1]**2
    f[1] = e*IonCharge/(u[1]**2)
    a[1] = f[1]/me
    
    for n in range(1,N):
        f[n+1] = e*IonCharge/(u[n]**2)
        a[n+1] = f[n+1]/me
        v[n+1] = v[n] +a[n+1]*t[n+1]
        u[n+1] = u[n] + v[n+1]*t[n+1] - .5*a[n+1]*t[n+1]**2
        
    #return u
    #return f
    return a
    #return v
    #print (u)
    #print(f)
    print(a)
    #print(v)
def vel_y(v_i, t, e, IonCharge, me, y_i):
    N = t.size - 1
    dt = t[1] - t[0]  
    u = np.zeros(N+1)
    v = np.zeros(N+1)
    f = np.zeros(N+1)
    a = np.zeros(N+1)
    u[0] = y_i
    f[0] = e*IonCharge/(u[0]**2)
    a[0] = f[0]/me
    v[0] = v_i
    v[1] = v[0] + a[1]*t[1]
    u[1] = u[0] + v[1]*t[1] - .5*a[1]*t[1]**2
    f[1] = e*IonCharge/(u[1]**2)
    a[1] = f[1]/me
    
    for n in range(1,N):
        f[n+1] = e*IonCharge/(u[n]**2)
        a[n+1] = f[n+1]/me
        v[n+1] = v[n] +a[n+1]*t[n+1]
        u[n+1] = u[n] + v[n+1]*t[n+1] - .5*a[n+1]*t[n+1]**2
        
    #return u
    #return f
    #return a
    return v
    #print (u)
    #print(f)
    #print(a)
    print(v)


y = position_y(0, time, e, IonCharge, me, yi_e)
print(y)
f_y = force_y(0, time, e, IonCharge, me, yi_e)
print(f_y)
a_y = acc_y(0, time, e, IonCharge, me, yi_e)
print(a_y)
v_y = vel_y(0, time, e, IonCharge, me, yi_e)
print(v_y)

plt.plot(time,y)
plt.xlabel("time")
plt.ylabel("y position")
plt.title("Y Position verus time")
plt.show()

#debug
plt.plot(time, f_y)
plt.xlabel("time")
plt.ylabel("y force")
plt.title("Y Force verus time")
plt.show()

plt.plot(time, v_y)
plt.xlabel("time")
plt.ylabel("y Velocity")
plt.title("Y Velocity verus time")
plt.show()

plt.plot(time,a_y)
plt.xlabel("time")
plt.ylabel("y acceleration")
plt.title("Y acceleration verus time")
plt.show()

def acc_r(v_ix, v_iy, t, e, IonCharge, me, y_i, x_i):
    N = t.size - 1
    dt = t[1] - t[0]  
    u = np.zeros(N+1)
    v = np.zeros(N+1)
    f = np.zeros(N+1)
    a = np.zeros(N+1)
    u[0] = math.sqrt(x_i**2 + y_i**2)
    f[0] = e*IonCharge/(u[0]**2)
    a[0] = f[0]/me
    v[0] = math.sqrt(v_ix**2+v_iy**2)
    v[1] = v[0] + a[1]*t[1]
    u[1] = u[0] + v[1]*t[1] - .5*a[1]*t[1]**2
    f[1] = e*IonCharge/(u[1]**2)
    a[1] = f[1]/me
    
    for n in range(1,N):
        f[n+1] = e*IonCharge/(u[n]**2)
        a[n+1] = f[n+1]/me
        v[n+1] = v[n] +a[n+1]*t[n+1]
        u[n+1] = u[n] + v[n+1]*t[n+1] - .5*a[n+1]*t[n+1]**2
        
    #return u
    #return f
    return a
    #return v
    #print (u)
    #print(f)
    print(a)
    #print(v)

radial_acc = acc_r(vi_e, 0, time, e, IonCharge, me, yi_e, xi_e)

plt.plot(time,radial_acc)
plt.xlabel("time")
plt.ylabel("r acceleration")
plt.title("Y acceleration verus time")
plt.show()

fourier = fft(radial_acc)

plt.plot(freq, fourier)
plt.ylabel("Power")
plt.xlabel("Frequency")
plt.title("Frequency verus Power Spectrum")
plt.show()

radial_acc1 = acc_r(vi_e, vi_e, time, e, IonCharge, me, yi_e, xi_e)

freq = [1/a for a in time]

fourier1 = fft(radial_acc1)

plt.plot(freq, fourier1)
plt.ylabel("Power")
plt.xlabel("Frequency")
plt.title("Frequency verus Power Spectrum, Vi_y = 10^7, b = 5*a0")
plt.show()

radial_acc2 = acc_r(vi_e, 0, time, e, IonCharge, me, 3*yi_e, xi_e)
fourier2 = fft(radial_acc2)
plt.plot(freq, fourier2)
plt.ylabel("Power")
plt.xlabel("Frequency")
plt.title("Frequency verus Power Spectrum, Vi_y = 0, b = 15*a0")
plt.show()

radial_acc3 = acc_r(vi_e, -vi_e, time, e, IonCharge, me, yi_e, xi_e)
fourier3 = fft(radial_acc3)
plt.plot(freq, fourier3)
plt.ylabel("Power")
plt.xlabel("Frequency")
plt.title("Frequency verus Power Spectrum, Vi_y = -10**7, b = 5*a0")
plt.show()

