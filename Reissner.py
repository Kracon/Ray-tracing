import numpy as np
import math
import scipy.integrate as sci
import matplotlib.pyplot as plt
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

def Reis(t, v): #Reissner-Nordström equations of motion.
    return np.array([v[2], L/v[0]**2, -L**2 * (3*M*v[0]-v[0]**2-2*Q)/v[0]**5]) 

def R(t, v): #The celestial sphere of the observer.
    return (v[0]*math.cos(v[1]) + r0)**2 + (v[0]*math.sin(v[1]))**2 - rlimit**2

R.terminal = True #Finishes at the celestial sphere.

def horizon(t, v): #The event horizon of the black hole.
    return v[0]-(M+np.sqrt(M**2-Q)) 

horizon.terminal = True #Finishes at event horizon.

def hole(radius, shift=[0,0]): #Plots the event horizon filled in black.
    C = [[shift[0]+radius*math.cos(x), shift[1]+radius*math.sin(x)] for x in np.linspace(0,2*math.pi,80)]
    X,Y = zip(*C)
    plt.plot(X,Y, 'black')
    plt.fill(X,Y, 'black')
    
def Celestial(radius, shift=[0,0]): #Plots the celestial sphere.
    C = [[shift[0]+radius*math.cos(x), shift[1]+radius*math.sin(x)] for x in np.linspace(0,2*math.pi,80)]
    X,Y = zip(*C)
    plt.plot(X,Y, 'black')
    
def orbitpath(r, phi): #Plots the solution of the equations of motion.
    polar = zip(r, phi)
    cart = map(lambda x: [x[0]*math.cos(x[1]), x[0]*math.sin(x[1])], polar)
    x,y = zip(*cart)
    p = plt.plot(x,y, 'r-')
    plt.axis([-3,3,-3,3])
    hole(M+np.sqrt(M**2-Q))
    plt.axes().set_aspect('equal')
    return p

def orbitpath2(r, phi): #Same as first but for different line colour/ types.
    polar = zip(r, phi)
    cart = map(lambda x: [x[0]*math.cos(x[1]), x[0]*math.sin(x[1])], polar)
    x,y = zip(*cart)
    p = plt.plot(x,y,'k--')
    plt.axis([-10,10,-10,10])
    hole(M+np.sqrt(M**2-Q))
    plt.axes().set_aspect('equal')
    return p

def f(x): #Used to interpolate apparent angle.
    return np.interp(x,alpha1,alpha1s)

def pixelmap(k,l): #Uses relationship to map pixels.
    i = k - blackhole[0]
    j = l - blackhole[1]
    r = np.sqrt((i)**2+(j)**2)
    theta = np.arctan2(j,i)
    phi = np.arctan2(r,fclength)
    x,y = (np.tan(f(phi))*np.cos(theta)*fclength, np.tan(f(phi))*np.sin(theta)*fclength)
    if -blackhole[0] <= x <= im.size[0]-blackhole[0] and -blackhole[1] <= y <= im.size[1]-blackhole[1]:
        return pixels[x+blackhole[0],y+blackhole[1]]
    else:
        return (0,0,0)

im = Image.open("Nebula_3.png")
new_im = Image.new('RGB', (im.size[0],im.size[1]), 'black')

trange = [0., 300.] #Time range.
fclength = im.size[0]/3 #Sets focal length of the camera.
blackhole = [im.size[0]/2,im.size[1]/2] #Sets location of black hole in the image.
M = 1 #Mass of the black hole.
Q = 9*M**2/8 #Charge of the black hole, using notation choice Q = {r_Q}^2.
r0 = 30 #Radius of where rays originate from.
rlimit = 200 #Distance we have set celestial sphere from oberser.
phi0 = np.pi #Angle from which rays originate from.
alpha1 = []
alpha1s = []
alpha2 = []
alpha2s = []
#Initial conditions.
v0p1 = [[[r0,math.pi,-1], L] for L in np.linspace(10.1, 100.0, 35)]
v0p2= [[[r0,math.pi,-1], L] for L in np.linspace(5.2,10, 20)]
v0p = v0p2 + v0p1

#v0s1 = [[[r0,0.45,-1], 14.6]]
#v0s2 = [[[r0,0.35,-1], 10.8]]
#v0s3 = [[[r0,0.25,-1], 7.6]]
#v0s4 = [[[r0,0.55,-1], 18.5]]
#v0p = v0s1 + v0s2 + v0s3 + v0s4

#v0p = [[[r,0,0],6] for r in np.linspace(4,17,5)]
#v0p = [[[10,np.pi,-1],5.8685]]
#v0p = [[[r,0,0],1] for r in np.linspace((3*M + np.sqrt(9*M**2 - 8*Q))/2 - 0.1,(3*M + np.sqrt(9*M**2 - 8*Q))/2 + 0.1,2)]

for v0 in v0p: #Solves equations of motion in range of initial conditions.
    L = v0[1]
    sol = sci.solve_ivp(Reis, trange, np.array(v0[0]), events = [horizon, R], t_eval = np.linspace(trange[0], trange[1], 10000))
    p = orbitpath(sol.y[0], sol.y[1])
    rfinal = sol.y[0][-1]
    phifinal = sol.y[1][-1]
    xfinal = rfinal * math.cos(phifinal)
    yfinal = rfinal * math.sin(phifinal)
    alpha1.append(math.atan2(v0[1]/v0[0][0], -v0[0][2])) #Appends apparent angle.
    alpha1s.append(-math.atan2(yfinal, r0 + xfinal)) #Appends actual angle.
Celestial(rlimit, shift = [-r0,0])
plt.xlabel('x/M')
plt.ylabel('y/M')
#plt.savefig('slightE.png', dpi=300) #Used to save images.
#plt.show()

#For comparison to Schwarzschild case.

#v0p1 = [[[r0,math.pi,-1], L] for L in np.linspace(10.1, 100.0, 35)]
#v0p2= [[[r0,math.pi,-1], L] for L in np.linspace(5.8,10, 20)]
#v0p = v0p2 + v0p1
#
for v0 in v0p: #Solves equations of motion in range of initial conditions.
    Q=0
    L = v0[1]
    sol = sci.solve_ivp(Reis, trange, np.array(v0[0]), events = [horizon, R], t_eval = np.linspace(trange[0], trange[1], 10000))
    p = orbitpath2(sol.y[0], sol.y[1])
    rfinal = sol.y[0][-1]
    phifinal = sol.y[1][-1]
    xfinal = rfinal * math.cos(phifinal)
    yfinal = rfinal * math.sin(phifinal)
    alpha2.append(math.atan2(v0[1]/v0[0][0], -v0[0][2])) 
    alpha2s.append(-math.atan2(yfinal, r0 + xfinal))
Celestial(rlimit, shift = [-r0,0])
plt.xlabel('x/M')
plt.ylabel('y/M')
plt.plot(0,0,'r-',label = 'Reissner-Nordström Black Hole')
plt.plot(0,0,'k--', label = 'Schwarzschild Black Hole')
plt.legend(loc = 'lower right')
plt.savefig('LoopR.png', dpi=300) #Used to save images.
plt.show()


#Produces the angle map.
n = np.linspace(0,1.4,100)
plt.plot(alpha1,alpha1s, label ='Reissner-Nordström Black Hole')
#plt.plot(alpha2,alpha2s, '--', label ='Schwarzschild Black Hole')
plt.plot(n,n,'k--')
plt.xlabel('Alpha Prime')
plt.ylabel('Alpha')
plt.axis([0,1.5,-1.5,1.5])
plt.legend(loc = 'lower right')
#plt.savefig('AngleMapR.png', dpi=300)
plt.show()

    
#pixels = im.load()
#new_pixels = new_im.load()
#
#for i in range(im.size[0]):
#    for j in range(im.size[1]):
#        new_pixels[i,j] = pixelmap(i,j)
#new_im.show()
#im.show()
#
#new_im.save("NebulaBlackHoleR.png", 'PNG') #Saves the produced image.