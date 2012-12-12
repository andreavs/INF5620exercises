from PDE_tools import *

def alpha1(u):
	a = 1.
	return a

def alpha2(u):
	a = 1+u**2
	return a

def alpha3(u):
	beta = 2.
	a = 1+beta*u*u
	return u



def run_main(divisions, degree, dt, T, alpha, rho, u0, u_e=None, f=None, savenumpy = False):
	# Physical constants:
	rho = 1.

	# Time discretisation:
	Nt = T/dt+1
	t = 0.0

	# FEniCS spesifics
	d = len(divisions)
	domain_type = [UnitInterval, UnitSquare, UnitCube]
	mesh = domain_type[d-1](*divisions)
	V = FunctionSpace(mesh, 'Lagrange', 1)
	#u0 = Expression('exp(-pi*pi*t)*cos(pi*x[0])',pi=pi,t=0)
	u0.t = 0
	u_p = project(u0, V)
	u = TrialFunction(V)
	v = TestFunction(V)
	if f is None:
		f = Constant(0.0)
	a = u*v*dx + dt/rho*alpha(u_p)*inner(nabla_grad(u), nabla_grad(v))*dx
	L = (u_p + dt/rho*f)*v*dx
	u = Function(V)
	A = assemble(a)
	u = Function(V)
	b = None
	c = None
	t = dt
	print savenumpy
	counter = 0
	while t<= T:
		u0.t = t
		f.t = t
		b = assemble(L,tensor=b)
		c = assemble(a,tensor=c)
		solve(a == L,u)
		u_p.assign(u)
		if u_e is not None: 
			u_e.t = t
			u_eplot = interpolate(u_e, V)
		filename = ''
		usave = np.zeros(1,'float')
		if savenumpy:
			usave, X,Y = numpyfy(u, mesh, divisions[0], divisions[1])
			filename = 'solution_%06d.txt' % counter
        	np.savetxt(filename, usave, fmt='%.18e', delimiter=' ', newline='\n')
		#u_e = Expression('t*x[0]*x[0]*(1./2 - x[0]/3.)', t = t)
		#u_etest = interpolate(u_e, V)
		#maxdiff = np.abs(u_etest.vector().array()-u.vector().array()).max()
		#print 'Max error, t=%.2f: %-10.16f' % (t,maxdiff)
		t += dt
		counter +=1
	if u_e is not None:
		return u_eplot, u
	else: 
		return u

def test_manufactured_solution():
	"""
	test for u = txx(1/2 - x/3), with alpha(u) = 1+uu, 
	"""
	rho = 1.0
	u0 = Constant('0.0')
	u_e = Expression('t*x[0]*x[0]*(1./2 - x[0]/3.)', t = 0)
	f = Expression('-rho*x[0]*x[0]*x[0]/3. + rho*x[0]*x[0]/2. + pow(t,3)*pow(x[0],4)*(pow(x[0],3)*8./9. - \
		28.*pow(x[0], 2)/9. + 7.*pow(x[0],1)/2. - 5./4.) + 2.*t*x[0] - t', t=0, rho=rho)
	
	dt = 0.1
	xnodes = 30
	T = [0.1, 0.5, 1.0, 2.0, 3.0]
	E = np.zeros(len(T), 'float')
	for i in range(len(T)):
		u_eplot, u = run_main([xnodes], 1, dt, T[i], alpha2, rho, u0, u_e, f)
		e = u_eplot.vector().array() - u.vector().array()
 		E[i] = np.sqrt(np.sum(e**2)/u.vector().array().size)
 		plt.plot(np.linspace(0,1,xnodes+1), np.abs(u.vector().array()-u_eplot.vector().array()))
 		#plt.plot(u.vector().array())
 		#plt.plot(u_eplot.vector().array())
 	return E

def test_convergence_rate():
	"""
	test the convergence rate in a spesific case. 
	"""
	h = 0.1
	T = 1.
	errorlist = []
	u0 = Expression('exp(-pi*pi*t)*cos(pi*x[0])',pi=pi,t=0)
	u_e = Expression('exp(-pi*pi*t)*cos(pi*x[0])',pi=pi,t=0)
	for i in range(7):
		dt = h
		xnodes = ynodes = int(1/sqrt(h))
		u_e, u = run_main([xnodes,ynodes], 1, dt, T, alpha1, 1.0, u0, u_e)
		e = u_e.vector().array() - u.vector().array()
 		E = np.sqrt(np.sum(e**2)/u.vector().array().size)
 		errorlist.append(E/h)
 		h = h/2
 	errordiff = []
 	for i in range(6):
 		errordiff.append(errorlist[i+1]/errorlist[i])
 	#	nt.assert_almost_equal(errordiff, 1, delta=0.1)
 	return errordiff

def test_convergence_rate2():
 	"""
 	test the convergence rate for a manufactured solution
 	"""
 	rho = 1.0
	u0 = Constant('0.0')
	u_e = Expression('t*x[0]*x[0]*(1./2 - x[0]/3.)', t = 0)
	h = 0.1
	T = 1.0
	errorlist = []
	for i in range(7):
		dt = h
		xnodes = ynodes = int(1/sqrt(h))
		f = Expression('rho*x[0]*x[0]*(-2.*x[0] + 3.)/6. - (-12.*t*x[0] + 3.*t*(-2.*x[0] + 3.))*(x[0]*x[0]*x[0]*x[0]*(-dt + t)*(-dt + t)*(-2.*x[0] \
			+ 3.)*(-2.*x[0] + 3.) + 36.)/324. - (-6.*t*x[0]*x[0] + 6.*t*x[0]*(-2.*x[0] + 3.))*(36.*pow(x[0],4)*(-dt + t)*(-dt + t)*(2.*x[0] - 3.) \
			+ 36.*x[0]*x[0]*x[0]*(-dt + t)*(-dt + t)*(-2.*x[0] + 3.)*(-2.*x[0] + 3.))/5832.', t=0, rho=rho, dt=dt)
		u_e, u = run_main([xnodes,ynodes], 1, dt, T, alpha1, 1.0, u0, u_e)
		e = u_e.vector().array() - u.vector().array()
 		E = np.sqrt(np.sum(e**2)/u.vector().array().size)
 		errorlist.append(E/h)
 		h = h/2
 	errordiff = []
 	for i in range(6):
 		errordiff.append(errorlist[i+1]/errorlist[i])
 	#	nt.assert_almost_equal(errordiff, 1, delta=0.1)
 	return errordiff

def gaussian():
	sigma = 0.5
	u0 = Expression('exp(-1./(2*sigma)*(x[0]*x[0]+ x[1]*x[1]))', sigma=sigma)
	alpha = alpha3
	dt = 0.01
	xnodes, ynodes = 10,10
	T = 2.0
	run_main([xnodes, ynodes], 1, dt,T, alpha, 1.0, u0, savenumpy=True)



if __name__ == '__main__':
	gaussian()
	mcrtmv(int(2./0.01-1), 0.01,1.0,1.0,11,11,savemovie=True, mvname='gaussian')
	# plt.hold('on')
	# print E
	# plt.legend(['Error for T = 0.1', 'Error for T = 0.5', 'Error for T = 1.0', 'Error for T = 3.0'], loc='upper left')
	# plt.title('Errors for different final times')
	# plt.ylabel('Absolute error')
	# plt.xlabel('position')
	# plt.show()