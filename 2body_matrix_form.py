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

		# x_and_v.append(particles.tolist())
		time += dt
	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0
	return e_error#x_and_v, e_error


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

		# x_and_v.append(particles.tolist())

		time += dt

	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return e_error #x_and_v, e_error

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

		# x_and_v.append(particles.tolist())

		time += dt
	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return e_error #x_and_v, e_error


def hermite(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0
	iterations = 2
	alpha = 7/6.

	while time < tfinal:
		print 'time =', time/pars.yr
		acc, jerk = fn.forces_hermite(particles, marr)

		old_x = np.copy(particles[:, 0, :])
		old_v = np.copy(particles[:, 1, :])
		old_a = np.copy(acc)
		old_j = np.copy(jerk)

		particles[:, 0, :] += particles[:, 1, :] * dt + acc * dt**2 / 2 + jerk * dt**3 / 6
		particles[:, 1, :] += acc * dt + jerk * dt**2 / 2
		for i in range(iterations):
			varr = old_v + (old_a + acc)*dt/2 + ((old_j - jerk)*dt**2)/12
			xarr = old_x + (old_v + particles[:, 1, :])*dt/2 + ((alpha/10)*((old_a - acc)*dt**2)) + ((6*alpha - 5)/120)*(old_j + jerk)*dt**3

		x_and_v.append(particles.tolist())
		print particles

		etot1 = fn.e_tot(particles, marr)
		e_error = (etot1 - etot0) / etot0
		return e_error, x_and_v


def plot(x1_val, y1_val, x2_val, y2_val):
	# plt.plot(x1_val, y1_val)
	# plt.plot(x2_val, y2_val)
	# print x2_val
	# plt.show()
	pass

def plottest(particles):
	'''
	Plots list (t, Np, 2, 3)
	'''
	xarr = np.array(particles)
	for planet in range(pars.Np):
		plt.plot(xarr[:, planet, 0, 0], xarr[:, planet, 0, 1])
	plt.show()

def plot_error(timestep, error1, error2, error3):
	plt.plot(timestep, np.abs(error1), label = 'Euler foward')
	plt.plot(timestep, np.abs(error2), label = 'midpoint')
	plt.plot(timestep, np.abs(error3), label = 'leapfrog')

	plt.legend()
	plt.yscale("log")
	plt.xscale("log")
	plt.show()


def main():
	error_hermite, pos_hermite = hermite(0.001*pars.yr, 3e8)
	# plottest(pos_hermite)
	# error_euler = []
	# error_midpoint = []
	# error_leapfrog =[]
	# timestep=[1e-2, 1e-3, 1e-4, 1e-5,1e-6]
	# for t in timestep:
	#
	# 	# print 'error euler forward    ' +str(timestep)+': ', euler   (timestep * pars.yr ,1e8)[1]
	# 	# print 'error midpoint forward ' +str(timestep)+': ', midpoint(timestep * pars.yr ,1e8)[1]
	# 	# print 'error leapfrog forward ' +str(timestep)+': ', leapfrog(timestep * pars.yr ,1e8)[1]
	# 	error_euler.append(euler   (t * pars.yr ,1e8).tolist())
	# 	error_midpoint.append(midpoint   (t * pars.yr ,1e8).tolist())
	# 	error_leapfrog.append(leapfrog   (t * pars.yr ,1e8).tolist())

	# print 'error_euler: ', error_euler
	# print 'error_midpoint: ', error_midpoint
	# print 'error_leapfrog: ', error_leapfrog
	# plot_error(timestep, error_euler, error_midpoint, error_leapfrog)


	# pos_leapfrog, error_leapfrog = leapfrog(dt, tfinal)
	# plottest(pos_leapfrog)
	# print 'error leapfrog =' ,error_leapfrog
	# print xv

if __name__ == '__main__':
	main()
