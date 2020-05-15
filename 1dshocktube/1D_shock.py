import numpy as np
import matplotlib.pyplot as plt
import csv


with open('shocktube_openfoam.csv', 'rU') as f:
    reader = csv.reader(f)
    header = next(reader)
    e = list(reader)
    f2= np.asarray(e) 

p5 = [float(x) for x in f2[:,0]]
x5 = [float(x) for x in f2[:,13]] 
r5 = [float(x) for x in f2[:,6]]
T5 = [float(x) for x in f2[:,4]]
u5 = [float(x) for x in f2[:,1]] 
c5 = [float(x) for x in f2[:,5]] 

p55 = []
x55 = []
r55 = []
T55 = []
u55 = []
c55 = []
mach = []
for i in range(len(p5)):
    p55.append(p5[i]/1e6)
    x55.append(x5[i])
    r55.append(r5[i])
    T55.append(T5[i])
    u55.append(u5[i])
    c55.append(c5[i])
    mach.append(np.sqrt(u5[i]*u5[i])/c5[i])
    
d = np.loadtxt("SIN32.0008")
e1 = np.loadtxt("1d_shocktube_press.tsv")
e2 = np.loadtxt("1d_shocktube_temp.tsv")
e3 = np.loadtxt("1d_shocktube_r.tsv")
e4 = np.loadtxt("1d_shocktube_velo.tsv")


p3=[]
x3=[]
u3=[]
r3=[]
y3=[]

p3=d[:,0]/1e6
u3=d[:,7]
r3=d[:,1]
x3=d[:,8]
y3=d[:,3]
T3=d[:,4]

#
p1=[]
x1=[]
p5=e1[:,1]
x5=e1[:,0]
#r5=e1[:,2]

T5=e2[:,1]
xT=e2[:,0]

r5=e3[:,1]
xr=e3[:,0]

u5=e4[:,1]
xu=e4[:,0]

    
plt.figure(figsize=(7,4)) 
plt.subplot(221)
plt.plot(x3,p3,linewidth = 1,linestyle='-',color = 'blue')
plt.plot(x55,p55,linewidth = 1,linestyle='-.',color = 'red')
plt.plot(x5,p5,linewidth = 1,linestyle=':',color = 'black') 
plt.ylabel('Pressure $[MPa]$',fontsize=9,position=(0,0.5),rotation = "vertical")
plt.xlabel('$x$ [m]',fontsize=12,position=(0.5,1),rotation = "horizontal")
plt.axis([0,100,0.9,3.1])
plt.grid(True)
#  
plt.subplot(222)
plt.plot(x3,r3,linewidth = 1,linestyle='-',color = 'blue')
plt.plot(x55,r55,linewidth = 1,linestyle='-.',color = 'red') 
plt.plot(xr,r5,linewidth = 1,linestyle=':',color = 'black') 
plt.ylabel('Density $[kg.m^{-3}]$',fontsize=9,position=(0,0.5),rotation = "vertical")
plt.xlabel('$x$ [m]',fontsize=12,position=(0.5,1),rotation = "horizontal")
plt.axis([0,100,17,70])
plt.grid(True)


#
plt.subplot(223)
plt.plot(x3,u3,linewidth = 1,linestyle='-',color = 'blue') 
plt.plot(x55,u55,linewidth = 1,linestyle='-.',color = 'red') 
plt.plot(xu,u5,linewidth = 1,linestyle=':',color = 'black') 
plt.ylabel('Velocity $[m.s^{-1}]$',fontsize=9,position=(0,0.5),rotation = "vertical")
plt.xlabel('$x$ [m]',fontsize=12,position=(0.5,1),rotation = "horizontal")
plt.axis([0,100,-2,120])
plt.grid(True)

###
plt.subplot(224)
plt.plot(x3,T3,linewidth = 1,linestyle='-',color = 'blue',label='Fang et al. (2018)')
plt.plot(xT,T5,linewidth = 1,linestyle=':',color = 'black',label='Giljarhus et al. (2012)') 
plt.plot(x55,T55,linewidth = 1,linestyle='-.',color = 'red',label='rhoCentralFoam') 
plt.ylabel('Temperature [K]',fontsize=9,position=(0,0.5),rotation = "vertical")
plt.xlabel('$x$ [m]',fontsize=12,position=(0.5,1),rotation = "horizontal")
plt.axis([0,100,250,350])
plt.grid(True)
plt.tight_layout()
ax = plt.gca()
ax.legend(loc='upper center', bbox_to_anchor=(-0.2, -0.34), ncol=5,frameon=False)
#
#####################################
plt.xticks(size = 12)
plt.yticks(size = 12)


plt.savefig("SINTEFshocktube.pdf")
plt.show()