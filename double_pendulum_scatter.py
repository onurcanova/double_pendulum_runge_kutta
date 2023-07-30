import math as mt
import matplotlib.pyplot as plt
#---initial conditions for double pendulum---
m1 = 80
m2 = 2.0
L1 = 3.0
L2 = 1.8
g = 9.81
h = 0.001
th10 = mt.pi / 1.8
th20 = mt.pi / 2.0
w10 = 0.
w20 = 0.
t0 = 0.
dth0 = th10 - th20
th1 = [th10]
th2 = [th20]
w1 = [w10]
w2 = [w20]
t = [t0]
dth = [dth0]
y10 = -L1 * mt.cos(th10)
y20 = -L1 * mt.cos(th10) - L2 * mt.cos(th20)
x10 =  L1 * mt.sin(th10)
x20 =  L1 * mt.sin(th10) + L2 * mt.sin(th20)
x1_list = [x10]
y1_list = [y10]
x2_list = [x20]
y2_list = [y20]


# energy func to obtain energy loss
def energy (th1,th2,w1,w2):
    y1 = -L1 * mt.cos(th1)
    y2 = -L1 * mt.cos(th1) - L2 * mt.cos(th2)
    x1 =  L1 * mt.sin(th1)
    x2 =  L1 * mt.sin(th1) + L2 * mt.sin(th2)
    r1 =   (x1 ** 2 + y1 ** 2) ** 0.5
    r2 =   (x2 ** 2 + y2 ** 2) ** 0.5
    pe = m1 * g * y1 + m2 * g * y2
    ke1 = 0.5*m1*(r1*w1)**2 + 0.5*m2*(r2*w2)**2    
    ke2= 0.5*m1*(L1**2)*(w1**2) + 0.5*m2*((L1**2)*(w1**2)+(L2**2)*(w2**2)+2*L1*L2*w1*w2*mt.cos(th1-th2))
    total_energy1 = pe + ke1
    total_energy2 = pe + ke2
    return [pe,ke1,ke2,total_energy1,total_energy2]
pe_list = [energy(th10,th20,w10,w20)[0]]
ke1_list = [energy(th10,th20,w10,w20)[1]]
total_energy_list1 = [energy(th10,th20,w10,w20)[3]]
ke2_list = [energy(th10,th20,w10,w20)[2]]
total_energy_list2 = [energy(th10,th20,w10,w20)[4]]
constant_energy = [energy(th10,th20,w10,w20)[0]]

def fth1(th1,th2,w1,w2,t):
    return w1

def fth2(th1,th2,w1,w2,t):
    return w2

def fw1(th1,th2,w1,w2,t):
    c1 = mt.cos(dth[-1])
    c2 = mt.sin(dth[-1])
    c3 = mt.sin(th1)
    c4 = mt.sin(th2)
    
    a = (m1 + m2) * L1
    b = m2 * L2 * c1
    c = m2 * L2 * c1
    d = m2 * L2
    e = -m2 * L2 * w2 * w2 * c2 - (m1 + m2) * g * c3
    f = m2 * L1 * w1 * w1 * c2 - m2 * g * c4
    return (e * d - b * f) / (a * d - b * c)

def fw2(th1,th2,w1,w2,t):
    c1 = mt.cos(dth[-1])
    c2 = mt.sin(dth[-1])
    c3 = mt.sin(th1)
    c4 = mt.sin(th2)

    a = (m1 + m2) * L1
    b = m2 * L2 * c1
    c = m2 * L2 * c1
    d = m2 * L2
    e = -m2 * L2 * w2 * w2 * c2 - (m1 + m2) * g * c3
    f = m2 * L1 * w1 * w1 * c2 - m2 * g * c4
    return (a * f - c * e) / (a * d - b * c)

# runge=kutta 4
for n in range(0,40000):
    k1_th1 = fth1(th1[-1],th2[-1],w1[-1],w2[-1],t[-1])
    k1_th2 = fth2(th1[-1],th2[-1],w1[-1],w2[-1],t[-1])
    k1_w1 = fw1(th1[-1],th2[-1],w1[-1],w2[-1],t[-1])
    k1_w2 = fw2(th1[-1],th2[-1],w1[-1],w2[-1],t[-1])

    k2_th1 = fth1(th1[-1] + k1_th1 * h / 2,th2[-1] + k1_th2 * h / 2,w1[-1] + k1_w1 * h / 2,w2[-1] + k1_w2 * h / -2,t[-1] + h / 2)
    k2_th2 = fth2(th1[-1] + k1_th1 * h / 2,th2[-1] + k1_th2 * h / 2,w1[-1] + k1_w1 * h / 2,w2[-1] + k1_w2 * h / -2,t[-1] + h / 2)
    k2_w1 = fw1(th1[-1] + k1_th1 * h / 2,th2[-1] + k1_th2 * h / 2,w1[-1] + k1_w1 * h / 2,w2[-1] + k1_w2 * h / -2,t[-1] + h / 2)
    k2_w2 = fw2(th1[-1] + k1_th1 * h / 2,th2[-1] + k1_th2 * h / 2,w1[-1] + k1_w1 * h / 2,w2[-1] + k1_w2 * h / -2,t[-1] + h / 2)

    k3_th1 = fth1(th1[-1] + k2_th1 * h / 2,th2[-1] + k2_th2 * h / 2,w1[-1] + k2_w1 * h / 2,w2[-1] + k2_w2 * h / -2,t[-1] + h / 2)
    k3_th2 = fth2(th1[-1] + k2_th1 * h / 2,th2[-1] + k2_th2 * h / 2,w1[-1] + k2_w1 * h / 2,w2[-1] + k2_w2 * h / -2,t[-1] + h / 2)
    k3_w1 = fw1(th1[-1] + k2_th1 * h / 2,th2[-1] + k2_th2 * h / 2,w1[-1] + k2_w1 * h / 2,w2[-1] + k2_w2 * h / -2,t[-1] + h / 2)
    k3_w2 = fw2(th1[-1] + k2_th1 * h / 2,th2[-1] + k2_th2 * h / 2,w1[-1] + k2_w1 * h / 2,w2[-1] + k2_w2 * h / -2,t[-1] + h / 2)

    k4_th1 = fth1(th1[-1] + k3_th1 * h,th2[-1] + k3_th2 * h,w1[-1] + k3_w1 * h,w2[-1] + k3_w2 * h,t[-1] + h)
    k4_th2 = fth2(th1[-1] + k3_th1 * h,th2[-1] + k3_th2 * h,w1[-1] + k3_w1 * h,w2[-1] + k3_w2 * h,t[-1] + h)
    k4_w1 = fw1(th1[-1] + k3_th1 * h,th2[-1] + k3_th2 * h,w1[-1] + k3_w1 * h,w2[-1] + k3_w2 * h,t[-1] + h)
    k4_w2 = fw2(th1[-1] + k3_th1 * h,th2[-1] + k3_th2 * h,w1[-1] + k3_w1 * h,w2[-1] + k3_w2 * h,t[-1] + h)


    th1_new = th1[-1] + h / 6. * (k1_th1 + 2. * k2_th1 + 2. * k3_th1 + k4_th1)
    th2_new = th2[-1] + h / 6. *(k1_th2 + 2. * k2_th2 + 2. * k3_th2 + k4_th2)
    w1_new = w1[-1] + h / 6. * (k1_w1 + 2. * k2_w1 + 2. * k3_w1 + k4_w1)
    w2_new = w2[-1] + h / 6. * (k1_w2 + 2. * k2_w2 + 2. * k3_w2 + k4_w2)
    t_new = t[-1] + h
    dth_new = th1_new - th2_new

    th1.append(th1_new)
    th2.append(th2_new)
    w1.append(w1_new)
    w2.append(w2_new)
    t.append(t_new)
    dth.append(dth_new)

    pe_list.append(energy(th1[-1],th2[-1],w1[-1],w2[-1])[0])
    

    ke1_list.append(energy(th1[-1],th2[-1],w1[-1],w2[-1])[1])
    total_energy_list1.append(energy(th1[-1],th2[-1],w1[-1],w2[-1])[3])

    ke2_list.append(energy(th1[-1],th2[-1],w1[-1],w2[-1])[2])
    total_energy_list2.append(energy(th1[-1],th2[-1],w1[-1],w2[-1])[4])

    constant_energy.append(constant_energy[-1])

    y1_list.append( -L1*mt.cos(th1[-1]))
    y2_list.append(-L1*mt.cos(th1[-1])-L2*mt.cos(th2[-1]))
    x1_list.append(L1*mt.sin(th1[-1]))
    x2_list.append(L1*mt.sin(th1[-1]) + L2*mt.sin(th2[-1]))

 
'''
--------scatter options
plt.plot(t,th1,"red",label="theta1")
plt.plot(t,th2,"blue",label="theta2")
plt.plot(t,w1,"red",label="theta1")
plt.plot(t,w2,"blue",label="theta2")
plt.plot(th1,th2,"red",label="theta1")
plt.plot(w1,w2,"blue",label="theta2")

plt.scatter(th1,w1,s=0.01,c="red",label="theta1")
plt.scatter(th2,w2,s=0.01,c="blue",label="theta1")


plt.plot(t,constant_energy,".b",label="constant energy")
plt.plot(t,total_energy_list1,".r",label="total energy 1")
plt.plot(t,total_energy_list2,".y",label="total energy 2")
plt.plot(t,pe_list,".g",label="pe")
plt.plot(t,ke2_list,".r",label="ke2")
'''

plt.scatter(x1_list,y1_list,s=0.01,c="red",label="theta1")
plt.scatter(x2_list,y2_list,s=0.01,c="blue",label="theta2")

plt.grid()
plt.legend()
plt.show()
    