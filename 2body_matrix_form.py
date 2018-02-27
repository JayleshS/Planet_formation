import numpy as np
import pars
import functions_matrix_form as fn
import matplotlib.pyplot as plt

# xarr, varr, marr = fn.init_2body(0)
# xarr, varr, marr = fn.init_2body(0)
# etot0 = fn.e_tot(xarr, varr, marr)


# etot0 = fn.e_tot(xarr, varr, marr)

x_and_v =[]


def euler(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0

	while time < tfinal:
		acc = fn.forces(particles, marr)
		particles[:, 0, :] += particles[:, 1, :]*dt
		particles[:, 1, :] += acc*dt

		x_and_v.append(particles.tolist())

		time += dt
	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0
	return x_and_v, e_error


def midpoint(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0

	while time < tfinal:
		acc = fn.forces(particles, marr)
		particles_mid = np.zeros((pars.Np, 2, 3))
		particles_mid[:, 0, :] = particles[:, 0, :] + particles[:, 1, :] * dt / 2
		particles_mid[:, 1, :] = particles[:, 1, :] + acc * dt / 2

		amid = fn.forces(particles_mid, marr)
		particles[:, 0, :] += particles_mid[:, 1, :] * dt
		particles[:, 1, :] += amid * dt

		x_and_v.append(particles.tolist())

		time += dt

	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return x_and_v, e_error

def leapfrog(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0
	# print dingems
	while time < tfinal:
		acc   = fn.forces(particles, marr)
		particles[:,1,:] += acc* dt/2
		particles[:,0,:] += particles[:,1,:]*dt
		acc   = fn.forces(particles, marr)
		particles[:,1,:] += acc* dt/2

		x_and_v.append(particles.tolist())

		time += dt
	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return x_and_v, e_error


def hermite(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0
	iterations = 2

	while time < tfinal:
		acc, jerk = fn.forces_hermite(particles, marr)

		old_x = np.copy(particles[:, 0, :])
		old_v = np.copy(particles[:, 1, :])
		old_a = np.copy(acc)
		old_j = np.copy(jerk)

		''''check to see if dividing gives problems!'''
		particles[:, 0, :] += particles[:, 1, :] * dt + acc  * dt**2 / 2 + jerk * dt**3 / 6
		particles[:, 1, :] += acc                * dt + jerk * dt**2 / 2

		for i in range(iterations):
			acc, jerk = fn.forces_hermite(particles, marr)
			particles[:, 1, :] = old_v + (old_a + acc               ) * dt/2 + ((old_j - jerk)*dt**2)/12
			particles[:, 0, :] = old_x + (old_v + particles[:, 1, :]) * dt/2 + ((old_a - acc )*dt**2)/12

		x_and_v.append(particles.tolist())

		time += dt

	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return x_and_v, e_error


def plot(x1_val, y1_val, x2_val, y2_val):
	# plt.plot(x1_val, y1_val)
	# plt.plot(x2_val, y2_val)
	# print x2_val
	# plt.show()
	pass

def plottest(particles, plotlabel="plotlabel"):
	'''
	Plots list (t, Np, 2, 3)
	'''
	xarr = np.array(particles)
	for planet in range(pars.Np):
		plt.plot(xarr[:, planet, 0, 0], xarr[:, planet, 0, 1],label=str(plotlabel))
	plt.legend()
	plt.show()

	plot_error(timestep, error_euler, error_midpoint, error_leapfrog, error_hermite)

def plot_error(timestep, error1, error2, error3, error4):
	plt.plot(timestep, np.abs(error1), label = 'euler_forward')
	plt.plot(timestep, np.abs(error2), label = 'midpoint')
	plt.plot(timestep, np.abs(error3), label = 'leapfrog')
	plt.plot(timestep, np.abs(error4), label = 'hermite')

	plt.title('timestep = ' + str(timestep))
	plt.legend()
	plt.yscale("log")
	plt.xscale("log")
	plt.show()


# print 'euler_forward: \n',euler(1,10)[0]
#
# print 'euler_forward: \n',euler(1,10)[0]



def main():
	pos_hermite, error_hermite = hermite(0.001*pars.yr, 3*pars.yr)
	print error_hermite
	plottest(pos_hermite)

	error_euler    = []
	error_midpoint = []
	error_leapfrog = []
	error_hermite  = []

	timestep=[1e-2, 1e-3, 1e-4]
	for t in timestep:
		# error_hermite, pos_hermite = hermite(t, 30*pars.yr)

		# print 'error euler forward    ' +str(t)+': ', euler   (t * pars.yr ,10)
		# print 'error midpoint forward ' +str(t)+': ', midpoint(t * pars.yr ,1e8)[1]
		# print 'error leapfrog forward ' +str(t)+': ', leapfrog(t * pars.yr ,1e8)[1]

		error_euler.append(euler       (t * pars.yr, 3*pars.yr )[1])
		error_midpoint.append(midpoint (t * pars.yr, 3*pars.yr )[1])
		error_leapfrog.append(leapfrog (t * pars.yr, 3*pars.yr )[1])
		error_hermite.append (hermite  (t * pars.yr, 3*pars.yr )[1])

	# print 'error_euler: ', error_euler
	# print 'error_midpoint: ', error_midpoint
	# print 'error_leapfrog: ', error_leapfrog
	plot_error(timestep, error_euler, error_midpoint, error_leapfrog, error_hermite)


	# pos_leapfrog, error_leapfrog = leapfrog(dt, tfinal)
	# plottest(pos_leapfrog)
	# print 'error leapfrog =' ,error_leapfrog
	# print xv

if __name__ == '__main__':
	main()
