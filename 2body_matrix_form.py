import numpy as np
import pars
import functions_matrix_form as fn
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2   # Create thicker lines around plots
plt.rcParams['lines.linewidth'] = 2  # Create thicker lines in plots



def euler(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0
	x_and_v =[]


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
	x_and_v =[]


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


def leapfrog(dt, tfinal, tau, drag=False, init_e=0.0):
	particles, marr = fn.init_2body(init_e)
	etot0 = fn.e_tot(particles, marr)

	eccentricity = []
	semi_major_axis = []
	x_and_v =[]
	vkep_list = []
	time_list = []

	time = 0

	while time < tfinal:
		if drag:
			acc, v_kep = fn.forces_total(particles, marr, tau)
		else:
			acc = fn.forces(particles, marr)
			v_kep=[]
		particles[:,1,:] += acc* dt/2
		particles[:,0,:] += particles[:,1,:]*dt

		if drag:
			acc, v_kep = fn.forces_total(particles, marr, tau)
		else:
			acc = fn.forces(particles, marr)
		particles[:,1,:] += acc* dt/2

		x_and_v.append(particles.tolist())
		# print 'next'
		ecc, a = fn.get_orbital_elements(particles, marr)

		eccentricity.append(np.sqrt(sum(ecc)**2))
		semi_major_axis.append(a)
		vkep_list.append(v_kep)
		time_list.append(time)

		time += dt

	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return x_and_v, e_error, semi_major_axis, eccentricity, vkep_list, time_list

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


def plot_pos(particles):
	'''
	Plots list (t, Np, 2, 3)
	'''
	xarr = np.array(particles)

	for planet in range(pars.Np):
		plt.plot(xarr[:, planet, 0, 0], xarr[:, planet, 0, 1], label=planet)
	plt.scatter(0,0, c="y")
	plt.scatter(1.,0, c="b")

	plt.legend()
	plt.xlabel('$x_{pos}$[AU]')
	plt.ylabel('$y_{pos}$[AU]')
	# plt.xlim(-1.5, 1.5)
	# plt.ylim(-1.5, 1.5)

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
	dt     = 0.001
	tfinal = 10

	calc_step = 10
	omega_k = (2*np.pi)

	tau_vals= [5e1]
	# tau_vals = np.geomspace(5e-2, 5e2, num=11)
	for tau in tau_vals:
		print 'calculating', tau
		pos_leapfrog, _, _,_, v_kep, time = leapfrog(dt, tfinal, tau, drag=True)

		pos_arr = np.array(pos_leapfrog)
		diff = (pos_arr[:, 1, 0, :] - pos_arr[:, 0, 0, :])**2
		rji = np.sqrt(np.sum(diff, axis=1))
		dr = (rji[:-1]-rji[1:])/dt
		# plt.plot(time[:-1], dr/v_kep[:-1], label='{:0.2e}'.format(tau))

		np.save("vratio_dt=" + str(dt)+"_tfinal=" + str(tfinal) + "_tau=" + str(tau),  [np.array(time[:-1]),dr/v_kep[:-1]])

		# v_ratio = np.load("vratio_dt=" + str(dt)+"_tfinal=" + str(tfinal) + "_tau=" + str(tau)+".npy")

		# print v_ratio

		# plt.plot(v_ratio[0], v_ratio[1], label='{:0.2e}'.format(tau) )
		# print len (v_kep)
		# print time

		# a_array = np.array(a_leapfrog)
		# delta_a = (a_array[1::calc_step] - a_array[0:-1:calc_step])/(calc_step*dt)
		# delta_vkep = vkep[0::calc_step]
		# plt.plot(time, rji)
		# plot_pos(pos_leapfrog)
		# plt.show()
		# plt.plot(time, a_leapfrog, label='{:0.2e}'.format(tau))
		# plt.xscale("Log")
        # plt.yscale("Log")


	# plt.legend()
	# plt.show()

if __name__ == '__main__':
	main()
