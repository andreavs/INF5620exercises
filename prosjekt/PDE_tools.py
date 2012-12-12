import numpy as np
from dolfin import * 
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import os, sys
#from enthought.mayavi import mlab as ml
from mayavi import mlab as ml
import pylab as pl

def numpyfy(u, mesh, xnodes, ynodes):
	"""
	Creates the good numpy jazz from FEniCS data, plotting ready!
	"""
	tempu = u.vector().array()
	tempmesh = mesh.coordinates()
	Nx = xnodes + 1
	Ny = ynodes + 1
	x = tempmesh[0:Nx,0]
	y = tempmesh[0:-1:Nx,1]
	X,Y = np.meshgrid(x,y)
	umesh = np.zeros((Ny,Nx), 'float')
	for i in range(Ny):
		umesh[i,:] = tempu[Nx*i:Nx*(i+1)]
	return umesh, X,Y


def mcrtmv(frames, dt,Lx,Ly,Nx,Ny,savemovie=False, mvname='test'):
	x = np.linspace(0,Lx,Nx);
	y = np.linspace(0,Lx,Nx);
	X,Y = np.meshgrid(x,y);
	size = 500,500
	
	fig = ml.figure(size= size, bgcolor=(1.,1.,1.));

	#fig.scene.anti_aliasing_frames=07

	#extent = [0,Nx-1,0,Ny-1,-30,30]
	
	ml.clf(figure=fig)
	u = np.loadtxt('solution_%06d.txt'%1);
	fname = '_tmp%07d.png' % 1
	s = ml.surf(x,y,u,figure=fig,vmin=-1,vmax=1)
	ml.axes(extent=[0,Lx,0,Ly,-2,2])
	ml.colorbar()
	ml.xlabel('x position')
	ml.ylabel('y position')
	ml.zlabel('wave amplitude')
	if savemovie == True:
		pl.ion()
		arr = ml.screenshot()
		img = pl.imshow(arr)
		pl.axis('off')
	
	for i in range(2,frames):

		u = np.loadtxt('solution_%06d.txt'%i);
		s.mlab_source.scalars = u
		fname = '_tmp%07d.png' % i
		if savemovie == True:
			arr = ml.screenshot()
			img.set_array(arr)
			pl.savefig(filename=fname)#,figure=fig)
			print 'Saving frame', fname
			pl.draw()

	fig.scene.disable_render = False
	os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s.mpg" % mvname);

