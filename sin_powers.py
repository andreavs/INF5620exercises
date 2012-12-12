from approx1D import *
import numpy as np 
import sympy as sm
import mpmath as mp

x = sm.Symbol('x')
f = sm.sin(x)
phi = [x**(2*j + 1) for j in range(5)]
Omega = [0,0]
for k in range(2,13):
	Omega[1] = k*np.pi/2
	u = least_squares(f,phi,Omega)
	comparison_plot(f,u,Omega, filename="tmp{0}.pdf".format(k))
	figure()

ucoeff = mp.taylor(mp.sin,0.0,9)
u = 0.0
for i in range(len(ucoeff)):
	u = u + ucoeff[i]*x**i

comparison_plot(f,u,Omega, filename="taylor.pdf")