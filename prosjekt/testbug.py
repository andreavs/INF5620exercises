from PDE_tools import *

def alpha1(u):
	a = 1.
	return a

def alpha2(u):
	a = 1+u^2
	return a



def run_main(xnodes, ynodes, dt, T, alpha, rho, u0, u_e=None):
	# Physical constants:
	rho = 1.
	if alpha == 1:
		alpha = alpha1
	elif alpha == 2:
		alpha = alpha2


	# Time discretisation:
	Nt = T/dt+1
	t = 0.0

	# FEniCS spesifics
	mesh = UnitSquare(xnodes, ynodes)
	V = FunctionSpace(mesh, 'Lagrange', 1)
	u0 = Expression('exp(-pi*pi*t)*cos(pi*x[0])',pi=pi,t=0)
	u0.t = 0
	u_p = project(u0, V)
	u = TrialFunction(V)
	v = TestFunction(V)
	f = Constant(0.0)
	a = u*v*dx + dt/rho*alpha(u_p)*inner(nabla_grad(u), nabla_grad(v))*dx
	L = (u_p + dt/rho*f)*v*dx
	u = Function(V)
	A = assemble(a)
	u = Function(V)
	b = None
	c = None
	t = dt
	while t<= T:
		b = assemble(L,tensor=b)
		c = assemble(a,tensor=c)
		u0.t = t
		solve(a == L,u)
		t += dt
		u_p.assign(u)
		if u_e is not None: 
			u_e = interpolate(u0, V)
		#maxdiff = np.abs(u_e.vector().array()-u.vector().array()).max()
		#print 'Max error, t=%.2f: %-10.16f' % (t,maxdiff)
		#print 'Max error, t=%.2f: %-10.16f' % (t,maxdiff)
	return u_e, u


def test_convergence_rate():
	h = 0.1
	T = 1.
	errorlist = []

	u0 = Expression('exp(-pi*pi*t)*cos(pi*x[0])',pi=pi,t=0)
	u_e = Expression('exp(-pi*pi*t)*cos(pi*x[0])',pi=pi,t=0)
	for i in range(8):
		dt = h
		xnodes = ynodes = int(1/sqrt(h))
		u_e, u = run_main(xnodes, ynodes, dt, T, alpha1, 1.0, u0, u_e)
		e = u_e.vector().array() - u.vector().array()
 		E = np.sqrt(np.sum(e**2)/u.vector().array().size)
 		errorlist.append(E/h)
 		h = h/2
 	errordiff = []
 	for i in range(7):
Solving linear variational problem.
 		errordiff.append(errorlist[i+1]/errorlist[i])
 	#	nt.assert_almost_equal(errordiff, 1, delta=0.1)
 	return errordiff



if __name__ == '__main__':
	errordiff = test_convergence_rate()
	print errordiff