import numpy as np
from matplotlib import pylab as plt
import struct
import matplotlib.gridspec as gridspec
from PIL import Image, ImageEnhance
from math import pi
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import matplotlib.ticker as ticker
from matplotlib.patches import Ellipse

#Script to generate a video from multiple bin files (phytoplankton) 
#If there are frames missing --->> ffmpeg -f image2 -i %*.png out.avi

Lx = 120.        #X box dimension
Ly = 60.        #Y box dimension
###########FLOW PARAMETERS
ni=2.
L1 = Lx
A1=.203 ############## A1 = U/(L1**(1/3))
Nk = 1
phy = 0#np.pi/4. 
##########################

Li = np.zeros(Nk)
Ai = np.zeros(Nk)
ki = np.zeros(Nk)
ei = np.zeros(Nk)
wi = np.zeros(Nk)

for i in range(0,Nk):
	Li[i] = L1*(ni**-i)
	Ai[i] = A1*((Li[i])**(1./3.))
	ki[i] = 2.*pi/Li[i]
	ei[i] = Li[i]/10.
	wi[i] = ki[i]*Ai[i]
	
#Dumping
B = 1.1
nd = 1
#########

mesh = 0.2
Y, X = np.mgrid[0:Ly:complex(0,Ly/mesh), 0:Lx:complex(0,Lx/mesh)]
U, V = np.mgrid[0:0:complex(0,Ly/mesh), 0:0:complex(0,Lx/mesh)]
dump = 0.5*(np.tanh(Y - 2.) - np.tanh(Y - 58.))


##############################################################################################

def rst(U,V):
	U, V = np.mgrid[0:0:complex(0,Ly/mesh), 0:0:complex(0,Lx/mesh)]
	return U,V
	
def flow(U,V,X,Y,t):
	for i in range(0,1):
		U -= Ai[i]*np.cos(ki[i]*(X-ei[i]*np.sin(wi[i]*t)))*np.sin(ki[i]*(Y+phy))
		V +=  Ai[i]*np.sin(ki[i]*(X-ei[i]*np.sin(wi[i]*t)))*np.cos(ki[i]*(Y+phy))
	for i in range(1,Nk):
		U -= Ai[i]*np.cos(ki[i]*(X-ei[i]*np.sin(wi[i]*t)))*np.sin(ki[i]*(Y-dump*ei[i]*np.sin(wi[i]*t+phy)))
		V +=  Ai[i]*np.sin(ki[i]*(X-ei[i]*np.sin(wi[i]*t)))*np.cos(ki[i]*(Y-dump*ei[i]*np.sin(wi[i]*t+phy)))
	return U,V	
			
	

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

interval = 0.5
tut = []
pup = []

#File Name
os.system('ls "filenamehere" > list.txt')
gap = 4.
pp = 0
depth = 60
Lxx = 120
meshsize = 0.2


cm = plt.get_cmap("jet")

with open('list.txt') as f:
    lines = f.read().splitlines()
   
ii = 0
for s in lines:
	fil = str(s)

	t = fil.split('T')
	t = t[1].split('_U')
	t = float(t[0])
	Uu = fil.split('U_')
	Uu = Uu[1].split('D')
	Uu = float(Uu[0])
	D = fil.split('D')
	D = float(D[1])
	Uu = Uu

###################################################

	if(t>500.):
		continue
	if(t % 2. != 0):
		continue

	U,V = rst(U,V)
	U,V = flow(U,V,X,Y,0.)


	fig, axs = plt.subplots(figsize=(8, 4))
	
	A = np.fromfile(fil, dtype=np.float32)
	B = []
	c = 0
	for lines in A:
		if c > 1:
			B.append(lines)
		c = c+1
	B = np.array(B)
	A = B.reshape([int(depth/meshsize), int(Lxx/meshsize)])

	
	dn = 2 ;  
	aspect = 20
	pad_fraction = 0.5


	divider = make_axes_locatable(axs)
	width = axes_size.AxesY(axs, aspect=1./aspect)
	pad = axes_size.Fraction(pad_fraction, width)
	cax = divider.append_axes("right", size=width, pad=pad)

	avg = np.mean(A)
	sd = np.std(A)
	contrfac = .9
	speed = np.sqrt(V**2 + U**2)
	ln = 1.5*speed / speed.max()
	p2 = axs.imshow(A/avg, cmap = 'jet', origin="lower", interpolation = 'bicubic', vmin = max(0, 1. - contrfac*sd/avg) , vmax = 1. + contrfac*sd/avg, extent=[0,Lx,0,Ly])
	axs.streamplot(X, Y, V, U, density = dn, linewidth=ln, color='k')
	p22 = axs.contour(X, Y, A/avg, [1], colors='white', linewidths=2.5)
#	axs.clabel(p22, fontsize=9, inline=True)
	axs.set_ylabel('z (m)')
	axs.set_xlabel('x (m)')
	clb2 = fig.colorbar(p2, ax = axs, cax=cax)
	clb2.set_label('$\Theta$ / <$\Theta$>', rotation=90)
	axs.set(xlim=(0, Lx), ylim=(Ly, 0))
	

	
	plt.savefig(str(int(10.*t/(10.*gap))).zfill(4) + '.png', bbox_inches = 'tight')
	plt.close()
	pp += 1
	print ("t = " + str(t) + ", img = " + str(int(10.*t/(10.*gap))).zfill(4) + '.png')


vidname = 'U' + str(Uu)[:4] + 'D' + str(D*100./36.)[:5] + '.mp4'
os.system('ffmpeg -r 8 -f image2 -i %*.png ' + vidname)
os.system('rm *.png')
os.system('rm list.txt')
