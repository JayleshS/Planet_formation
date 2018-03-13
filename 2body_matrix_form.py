import numpy as np
import pars
import functions_matrix_form as fn
import matplotlib.pyplot as plt
# import matplotlib as mpl

plt.rcParams['axes.linewidth'] = 2   # Create thicker lines around plots
plt.rcParams['lines.linewidth'] = 2  # Create thicker lines in plots

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


def leapfrog(dt, tfinal, drag=False):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0
	# print dingems
	while time < tfinal:
		# print particles
		if drag:
			acc = fn.forces_total(particles, marr)
		else:
			acc = fn.forces(particles, marr)
		particles[:,1,:] += acc* dt/2
		particles[:,0,:] += particles[:,1,:]*dt

		if drag:
			acc = fn.forces_total(particles, marr)
		else:
			acc = fn.forces(particles, marr)
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
	ecc = []
	a_list = []

	while time < tfinal:
		acc, jerk = fn.forces_hermite(particles, marr)

		old_x = np.copy(particles[:, 0, :])
		old_v = np.copy(particles[:, 1, :])
		old_a = np.copy(acc)
		old_j = np.copy(jerk)

		particles[:, 0, :] += particles[:, 1, :] * dt + acc  * dt**2 / 2 + jerk * dt**3 / 6
		particles[:, 1, :] += acc                * dt + jerk * dt**2 / 2

		for i in range(iterations):
			acc, jerk = fn.forces_hermite(particles, marr)
			particles[:, 1, :] = old_v + (old_a + acc               ) * dt/2 + ((old_j - jerk)*dt**2)/12
			particles[:, 0, :] = old_x + (old_v + particles[:, 1, :]) * dt/2 + ((old_a - acc )*dt**2)/12

		x_and_v.append(particles.tolist())
		e, a = fn.get_orbital_elements(particles, marr)
		ecc.append(np.sqrt(sum(e**2)))
		a_list.append(a)

		time += dt

	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return x_and_v, e_error, ecc, a_list


def leapfrog_drag(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0

	while time < tfinal:
		

		acc   = fn.forces(particles, marr)[0]
		particles[0,1,:] += acc* dt/2
		particles[0,0,:] += particles[0,1,:]*dt
		acc   = fn.forces(particles, marr)[0]
		particles[0,1,:] += acc* dt/2

		grav   = fn.forces(particles, marr)[1]
		drag = fn.forces_migration(particles, marr)
		acc = grav + drag*1e-1

		particles[1,1,:] += acc* dt/2
		particles[1,0,:] += particles[1,1,:]*dt

		grav   = fn.forces(particles, marr)[1]
		drag = fn.forces_migration(particles, marr)
		acc = grav + drag*1e-1

		particles[1,1,:] += acc* dt/2

		x_and_v.append(particles.tolist())

		time += dt
	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return x_and_v, e_error


def plot_pos(particles):
	'''
	Plots list (t, Np, 2, 3)
	'''
	xarr = np.array(particles)
	for planet in range(pars.Np):
		plt.plot(xarr[:, planet, 0, 0], xarr[:, planet, 0, 1], label=planet)
	plt.axhline(y=pars.au)
	plt.legend()
	plt.xlabel('$x_{pos}$[cm]')
	plt.ylabel('$y_{pos}$[cm]')
	plt.show()


def plot_error(timestep, error1, error2, error3, error4):
	plt.plot(timestep, np.abs(error1), label = 'euler_forward')
	plt.plot(timestep, np.abs(error2), label = 'midpoint')
	plt.plot(timestep, np.abs(error3), label = 'leapfrog')
	plt.plot(timestep, np.abs(error4), label = 'hermite')

	plt.title('timestep = ' + str(timestep))
	plt.legend()
	plt.xlabel("stepsize [yr]")
	plt.ylabel("Error in energy")
	plt.yscale("log")
	plt.xscale("log")
	plt.show()


def main():
	pos_leapfrog, error_leapfrog = leapfrog(0.00001*pars.yr, 8*pars.yr, drag=True)
	plot_pos(pos_leapfrog)
	# xarr = np.array(pos_leapfrog)
	# plt.plot(xarr[:,1,0,0], xarr[:,1,0,1], label='earth')
	# plt.plot(xarr[:,0,0,0], xarr[:,0,0,1], label='sun')
	# plt.legend()
	# plt.show()
	# print 'error leapfrog =' ,error_leapfrog

if __name__ == '__main__':
	main()
