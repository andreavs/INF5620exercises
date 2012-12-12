from nonlinear_diffusion import *

def test_convergence_rate():
	h = 0.1
	T = 1.
	errorlist = []
	for i in range(4):
		dt = h
		xnodes = ynodes = int(1/sqrt(h))
		u_e = None
		u = None
		u_e, u = run_main(xnodes, ynodes, dt, T)
		e = u_e.vector().array() - u.vector().array()
 		E = np.sqrt(np.sum(e**2)/u.vector().array().size)
 		errorlist.append(E/h)
 		h = h/2

 	for i in range(3):
 		errordiff = errorlist[i+1]/errorlist[i]
 		nt.assert_almost_equal(errordiff, 1, delta=0.1)