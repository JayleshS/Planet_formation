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


def leapfrog(dt, tfinal, t_stop, drag=False, init_e=0.0):
	particles, marr = fn.init_2body(init_e)
	etot0 = fn.e_tot(particles, marr)

	eccentricity = []
	semi_major_axis = []
	x_and_v =[]
	vkeps = []
	time_list = []

	time = 0

	while time < tfinal:
		if drag:
			acc, vKep = fn.forces_total(particles, marr, t_stop)
		else:
			acc = fn.forces(particles, marr)
			vKep=[]
		particles[:,1,:] += acc* dt/2
		particles[:,0,:] += particles[:,1,:]*dt

		if drag:
			acc, vKep = fn.forces_total(particles, marr, t_stop)
		else:
			acc = fn.forces(particles, marr)
		particles[:,1,:] += acc* dt/2

		x_and_v.append(particles.tolist())
		# ecc, a = fn.get_orbital_elements(particles, marr)
		# if a < 0.046*pars.au:
		# 	break

		eccentricity.append(np.sqrt(sum(ecc)**2))
		semi_major_axis.append(a)
		vkeps.append(vKep)
		time_list.append(time)

		time += dt

	# print a/pars.au

	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return x_and_v, e_error, semi_major_axis, eccentricity, vkeps, time_list

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
	# tstop = 1*pars.yr
	dt = 0.001*pars.yr
	tfinal = 10*pars.yr    # v_head = vKep *0.4

	calc_step = 10
	omega_k = (2*np.pi)/pars.yr


	t_stop_factors= [5e-6]#,1e-6,1e-7,1e-8]
	# t_stop_factors=[2*np.pi/pars.yr]
	# t_stop_factors = np.geomspace(1e-5, 1e-8, num=8)
	for tstop in t_stop_factors:
		print 'calculating', tstop

		# plt.title("tstop_factor = "+ str(tstop))
		pos_leapfrog,_,a_leapfrog,_,_,time = leapfrog(dt, tfinal, tstop, drag=True)

		# pos_leapfrog, error_leapfrog, a_leapfrog, e_leapfrog, vkep = leapfrog(dt, tfinal, tstop, drag=True)

		# a_array = np.array(a_leapfrog)
		# delta_a = (a_array[1::calc_step] - a_array[0:-1:calc_step])/(calc_step*dt)
		# delta_vkep = vkep[0::calc_step]
		# plt.plot(delta_a/delta_vkep)
		plot_pos(pos_leapfrog)
		# plt.show()
		# plt.plot(time/pars.yr, a_leapfrog, label='{:0.2e}'.format(tstop))
		# plt.xscale("Log")
  #       plt.yscale("Log")
		# plt.plot(vkep)
		# plt.plot(a_leapfrog, label=tstop)
		# plt.xscale("Log")
        # plt.yscale("Log")
	plt.legend()
	plt.show()

if __name__ == '__main__':
	main()
