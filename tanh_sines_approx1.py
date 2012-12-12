from approx1D import *
import subprocess                 # For issuing commands to the OS.
import os
import numpy as np 
import sympy as sm
import mpmath as mp

x = sm.Symbol('x')
s = 20.
f = sm.tanh(s*(x-np.pi))
Omega = [0,2*np.pi]

u = 0.0
for k in range(15):
	phi = [sm.sin((2*k+1)*x)]
	u = u + least_squares(f,phi,Omega)
	comparison_plot(f,u,Omega, filename="tmp{0}.png".format(str(k).zfill(3)))
	clf()



# Now that we have graphed images of the dataset, we will stitch them
# together using Mencoder to create a movie.  Each image will become
# a single frame in the movie.
#
# We want to use Python to make what would normally be a command line
# call to Mencoder.  Specifically, the command line call we want to
# emulate is (without the initial '#'):
# mencoder mf://*.png -mf type=png:w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o output.avi
# See the MPlayer and Mencoder documentation for details.
#

command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=10',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'output.avi')

#os.spawnvp(os.P_WAIT, 'mencoder', command)

print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
subprocess.check_call(command)

print "\n\n The movie was written to 'output.avi'"

print "\n\n You may want to delete *.png now.\n\n"




