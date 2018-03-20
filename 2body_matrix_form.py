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


def plot_pos(particles, ax_range=1.2, tfinal=None):
	'''
	Plots list (t, Np, 2, 3)
	'''
	xarr = np.array(particles)

	plt.plot(xarr[:, 0, 0, 0], xarr[:, 0, 0, 1], c='y')
	plt.plot(xarr[:, 1, 0, 0], xarr[:, 1, 0, 1], label='Test particle', c='r')

	plt.scatter(1.,0, c="b", label='Initial postition')
	plt.scatter(0,0, c="y", label='Sun')

	plt.legend()
	plt.xlabel('$x_{pos}$(au)')
	plt.ylabel('$y_{pos}$(au)')
	plt.xlim(-ax_range, ax_range)
	plt.ylim(-ax_range, ax_range)
	plt.title('Total integration time: %s years' %tfinal)
	plt.axis('equal')
	plt.show()


def plot_element(time, element, xlabel='Time', ylabel=None):
	pass


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
		a_leapfrog = np.array(a_leapfrog)
		diff = (pos_arr[:, 1, 0, :] - pos_arr[:, 0, 0, :])**2
		rji = np.sqrt(np.sum(diff, axis=1))
		vr = (rji[:-1]-rji[1:])/dt

		# vr_a = (a_leapfrog[:-1] - a_leapfrog[1:])/dt
		# plt.plot(time[:-1], vr_a/v_kep[:-1])
		# plt.plot(time[:-1], vr/v_kep[:-1], label='{:0.2e}'.format(tau))
		#
		# np.save("vratio_dt=" + str(dt)+"_tfinal=" + str(tfinal) + "_tau=" + str(tau),  dr/v_kep[:-1])
		# print len(time[:-1])
		# print len(vr/v_kep[:-1])

		# ding = np.load("vratio_dt=0.001_tfinal=100_tau=1.0.npy")
		plot_pos(pos_leapfrog, tfinal=tfinal)
		# plt.plot(time, a_leapfrog, label='{:0.4e}'.format(tau))

		# print v_ratio

		# plt.plot(v_ratio[0], v_ratio[1], label='{:0.2e}'.format(tau) )
		# print len (v_kep)
		# print time

	# plt.legend()
	# plt.show()

if __name__ == '__main__':
	main()
