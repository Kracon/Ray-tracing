import numpy as np
import math
import scipy.integrate as sci
import matplotlib.pyplot as plt

def Kerr(t, v):
    return np.array([v[1],-(2*M*(v[2]**2 + a**2)*(v[2]**2 - a**2 *np.cos(v[4])**2))/(Sigma(v[2],v[4])**2*Delta(v[2]))*v[1]*v[3]+(4*M*a**2*v[2]*np.sin(v[4])*np.cos(v[4]))/(Sigma(v[2],v[4])**2)*v[1]*v[5]-(2*M*a*np.sin(v[4])**2 *(a**2 *np.cos(v[4])**2 *(a**2 - v[2]**2) -v[2]**2*(a**2 + 3*v[2]**2)))/(Sigma(v[2],v[4])**2 *Delta(v[2])) *v[3]*v[7]-(4*M*a**3 *v[2]*np.sin(v[4])**3 *np.cos(v[4]))/(Sigma(v[2],v[4])**2)*v[5]*v[7],v[3],-(M*Delta(v[2])*(v[2]**2 - a**2*np.cos(v[4])**2))/(Sigma(v[2],v[4])**3)*v[1]**2-(v[2]*a**2*np.sin(v[4])**2 - M*(v[2]**2 - a**2 *np.cos(v[4])**2))/(Sigma(v[2],v[4])*Delta(v[2]))*v[3]**2 +(2*Delta(v[2])*M*a*np.sin(v[4])**2*(v[2]**2 - a**2 *np.cos(v[4])**2))/(Sigma(v[2],v[4])**3)*v[1]*v[7]+(2*a**2*np.sin(v[4])*np.cos(v[4]))/(Sigma(v[2],v[4]))*v[3]*v[5]+(v[2]*Delta(v[2]))/(Sigma(v[2],v[4]))*v[5]**2-(Delta(v[2])*np.sin(v[4])**2*(M*a**2*np.sin(v[4])**2*(v[2]**2-a**2*np.cos(v[4])**2)-v[2]*Sigma(v[2],v[4])**2))/(Sigma(v[2],v[4])**3)*v[7]**2,v[5],(2*M*a**2*v[2]*np.sin(v[4])*np.cos(v[4]))/(Sigma(v[2],v[4])**3)*v[1]**2-(4*M*a*v[2]*(v[2]**2 + a**2)*np.sin(v[4])*np.cos(v[4]))/(Sigma(v[2],v[4])**3)*v[1]*v[7]-(a**2*np.sin(v[4])*np.cos(v[4]))/(Sigma(v[2],v[4])*Delta(v[2]))*v[3]**2-(2*v[2])/(Sigma(v[2],v[4]))*v[3]*v[5]+(a**2*np.sin(v[4])*np.cos(v[4]))/(Sigma(v[2],v[4]))*v[5]**2 +(np.sin(v[4])*np.cos(v[4])*((v[2]**2 + a**2)*Sigma(v[2],v[4])**2 + 2*M*a**2*v[2]*np.sin(v[4])**2*(Sigma(v[2],v[4]) + v[2]**2 + a**2)))/(Sigma(v[2],v[4])**3)*v[7]**2,v[7],-(2*M*a*(v[2]**2 - a**2 *np.cos(v[4])**2))/(Sigma(v[2],v[4])**2 * Delta(v[2]))*v[1]*v[3]+(4*M*a*v[2])/(np.tan(v[4])*Sigma(v[2],v[4])**2)*v[1]*v[5]-(2*(Sigma(v[2],v[4])**2 + 2*M*a**2*v[2]*np.sin(v[4])**2))/(np.tan(v[4])*Sigma(v[2],v[4])**2)*v[5]*v[7]-(2*v[2]*Sigma(v[2],v[4])**2+2*M*(a**4*np.sin(v[4])**2 *np.cos(v[4])**2 - v[2]**2 *(Sigma(v[2],v[4]) + v[2]**2 + a**2)))/(Sigma(v[2],v[4])**2*Delta(v[2]))*v[3]*v[7]])
 
def R(t, v): #The celestial sphere of the observer.
    return (v[2]*math.sin(v[4])*math.cos(v[6]) + r0)**2 + (v[2]*math.sin(v[4])*math.sin(v[6]))**2 + v[2]**2 * math.cos(v[4])**2 - rlimit**2

R.terminal = True #Finishes at the celestial sphere.

def horizon(t, v): #The event horizon of the black hole.
    return v[2]-(M+np.sqrt(M**2-a**2)) 

horizon.terminal = True #Finishes at event horizon.

def Sigma(r,theta):
    return r**2 + a**2 *np.cos(theta)**2

def Delta(r):
    return r**2 -2*M*r + a**2

def hole(radius, shift=[0,0]): #Plots the event horizon filled in black.
    C = [[shift[0]+radius*math.cos(x), shift[1]+radius*math.sin(x)] for x in np.linspace(0,2*math.pi,80)]
    X,Y = zip(*C)
    plt.plot(X,Y, 'black')
    plt.fill(X,Y, 'black')
    
def Celestial(radius, shift=[0,0]): #Plots the celestial sphere.
    C = [[shift[0]+radius*math.cos(x), shift[1]+radius*math.sin(x)] for x in np.linspace(0,2*math.pi,80)]
    X,Y = zip(*C)
    plt.plot(X,Y, 'black')
    
def orbitpath(r, phi, theta): #Plots the solution of the equations of motion.
    polar = zip(r, phi, theta)
    cart = map(lambda x: [np.sqrt(x[0]**2 + a**2)*math.cos(x[1])*math.sin(x[2]), np.sqrt(x[0]**2 + a**2)*math.sin(x[1])*math.sin(x[2])], polar)
    x,y = zip(*cart)
    p = plt.plot(x,y, 'r-') #, 'r-'
    plt.axis([-250,250,-250,250])
    hole(M+np.sqrt(M**2 - a**2))
    plt.axes().set_aspect('equal')
    return p

def orbitpath2(r, phi, theta): #Same as first but for different line colour/ types.
    polar = zip(r, phi, theta)
    cart = map(lambda x: [np.sqrt(x[0]**2 + a**2)*math.cos(x[1])*math.sin(x[2]), np.sqrt(x[0]**2 + a**2)*math.sin(x[1])*math.sin(x[2])], polar)
    x,y = zip(*cart)
    p = plt.plot(x,y, 'k--') #, 'r-'
    plt.axis([-250,250,-250,250])
    hole(M+np.sqrt(M**2 - a**2))
    plt.axes().set_aspect('equal')
    return p

def f(x): #Used to interpolate apparent angle.
    return np.interp(x,alpha1,alpha1s)

trange = [0., 300.] #Time range.
M = 1 #Mass of the black hole.
a = 0.99 #Angular momentum of the black hole.
rlimit = 200 #Distance we have set celestial sphere from oberser.
#Initial conditions
r0 = 30
rd0 = -1
th0 = np.pi/2
thd0 = 0
p0 = np.pi
pd0 = 0
alpha1 = []
alpha1s = []
#Initial conditions.
#v0s1 = [[0.45, 0.0163]]
#v0s2 = [[0.35, 0.012]]
#v0s3 = [[0.25, 0.0088]]
#v0s4 = [[0.55, 0.0204]]
#
#
#v0s5 = [[-0.45, -0.0163]]
#v0s6 = [[-0.35, -0.012]]
#v0s7 = [[-0.25, -0.0088]]
#v0s8 = [[-0.55, -0.0204]]
#
#pdd0 = v0s1 + v0s2 + v0s3 + v0s4 + v0s5 + v0s6 + v0s7 + v0s8
#
pd01 = [[L] for L in np.linspace(0.0045, 0.0075, 5)]
pd02 = [[L] for L in np.linspace(0.0075, 0.015, 7)]
pd03 = [[L] for L in np.linspace(0.015, 0.3, 30)]
pdd0 = pd01 + pd02 + pd03

#pd01 = [[L] for L in np.linspace(-pd0max, -0.015, 30)]
#pd02 = [[L] for L in np.linspace(-0.015, -0.008, 7)]
#pdd0 = pd01 + pd02
#pdd0 = [[L] for L in np.linspace(0.0065, 0.03, 4)]

for pd0s in pdd0: #Solves equations of motion in range of initial conditions.
    pd0 = pd0s[0]
    tdK0 = -(2*M*r0*a*np.sin(th0)**2)/(Sigma(r0, th0)-2*M*r0)*pd0 + np.sqrt((M**2 *r0**2 *a**2 *np.sin(th0)**4)/((Sigma(r0, th0)-2*M*r0)**2) *pd0**2 + Sigma(r0, th0)/(Sigma(r0, th0)-2*M*r0) *(Sigma(r0, th0)/Delta(r0)*rd0**2 + Sigma(r0, th0)*thd0**2 + (r0**2 + a**2 + 2*M*r0*a**2 /Sigma(r0, th0) * np.sin(th0)**2)*np.sin(th0)**2 *pd0**2))
    v0t = [0,tdK0,r0,rd0,th0,thd0,p0,pd0]
    sol = sci.solve_ivp(Kerr, trange, np.array(v0t), events = [horizon, R], t_eval = np.linspace(trange[0], trange[1], 10000))
    p = orbitpath(sol.y[2], sol.y[6], sol.y[4])
    rfinal = sol.y[2][-1]
    thetafinal = sol.y[4][-1]
    phifinal = sol.y[6][-1]
    xfinal = np.sqrt(rfinal**2 + a**2) * math.sin(thetafinal) * math.cos(phifinal)
    yfinal = np.sqrt(rfinal**2 + a**2) * math.sin(thetafinal) * math.sin(phifinal)
    zfinal = rfinal*math.cos(thetafinal)
    if yfinal >= 0: #Appends apparent angle.
        alpha1s.append(-math.atan2(np.sqrt(zfinal**2 + yfinal**2), r0 + xfinal))
    else:
        alpha1s.append(math.atan2(np.sqrt(zfinal**2 + yfinal**2), r0 + xfinal))

M = 0
a = 0

for pd0s in pdd0: #Solves equations of motion in range of initial conditions.
    pd0 = pd0s[0]
    tdK0 = -(2*M*r0*a*np.sin(th0)**2)/(Sigma(r0, th0)-2*M*r0)*pd0 + np.sqrt((M**2 *r0**2 *a**2 *np.sin(th0)**4)/((Sigma(r0, th0)-2*M*r0)**2) *pd0**2 + Sigma(r0, th0)/(Sigma(r0, th0)-2*M*r0) *(Sigma(r0, th0)/Delta(r0)*rd0**2 + Sigma(r0, th0)*thd0**2 + (r0**2 + a**2 + 2*M*r0*a**2 /Sigma(r0, th0) * np.sin(th0)**2)*np.sin(th0)**2 *pd0**2))
    v0t = [0,tdK0,r0,rd0,th0,thd0,p0,pd0]
    sol = sci.solve_ivp(Kerr, trange, np.array(v0t), events = [horizon, R], t_eval = np.linspace(trange[0], trange[1], 10000))
    p = orbitpath(sol.y[2], sol.y[6], sol.y[4])
    rfinal = sol.y[2][-1]
    thetafinal = sol.y[4][-1]
    phifinal = sol.y[6][-1]
    xfinal = np.sqrt(rfinal**2 + a**2) * math.sin(thetafinal) * math.cos(phifinal)
    yfinal = np.sqrt(rfinal**2 + a**2) * math.sin(thetafinal) * math.sin(phifinal)
    zfinal = rfinal*math.cos(thetafinal)
    alpha1.append(-math.atan2(yfinal, r0 + xfinal))
Celestial(rlimit, shift = [-r0,0])
plt.xlabel('x/M')
plt.ylabel('y/M')
#plt.savefig('slightE.png', dpi=300) #Used to save images.
plt.show()


#Produces the angle map.
n = np.linspace(0,1.4,100)
plt.plot(alpha1,alpha1s)
plt.plot(n,n,'k--')
plt.xlabel('Alpha Prime')
plt.ylabel('Alpha')
plt.axis([0,1.5,-1.5,1.5])
plt.legend(loc = 'lower right')
#plt.savefig('AngleMapK.png', dpi=300)
plt.show()