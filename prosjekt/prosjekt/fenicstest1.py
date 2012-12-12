import numpy as np
from dolfin import * 
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm

mesh = UnitSquare(6,4)
V = FunctionSpace(mesh, 'Lagrange', 1)
u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')

def u0_boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, u0, u0_boundary)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(nabla_grad(u), nabla_grad(v))*dx
L = f*v*dx


u = Function(V)
solve(a == L, u, bc)

#plot(u)
#plot(mesh)

#file = File('poisson.pvd')
#file << u.vector().array() 
# test = u.vector().array()
# printmesh = mesh.coordinates()
# print test
# x = printmesh[0:7,0]
# y = printmesh[0:-1:7,1]
# X, Y = np.meshgrid(x, y)
# meshtest = np.zeros((7,5), 'float')
# for i in range(5):
# 	meshtest[:,i] = test[7*i:7*i+7]

#a = u_box.grid.coorv[X]
#print printmesh


def numpyfy(u, mesh, xnodes, ynodes):
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

umesh, X, Y = numpyfy(u,mesh,6,4)

print umesh
print X

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
plt.figure()
CS = plt.contour(X,Y, umesh)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('Simplest default with labels')

plt.show()