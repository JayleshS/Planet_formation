import numpy as np
import pars
import functions_matrix_form as fn
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.optimize import curve_fit


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
	# save_file = open("object_time_x,v,m_dt=" + str(dt) + "_tfinal=" + str(tfinal) + "_tau=" +str(tau)+".csv", "w")

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
			acc, v_kep, v_r = fn.forces_total(particles, marr, tau)
		else:
			acc = fn.forces(particles, marr)
			v_kep=[]

		particles[:,1,:] += acc* dt/2
		particles[:,0,:] += particles[:,1,:]*dt

		if drag:
			acc, v_kep, v_r = fn.forces_total(particles, marr, tau, v_r=v_r)
		else:
			acc = fn.forces(particles, marr)
		particles[:,1,:] += acc* dt/2

		x_and_v.append(particles.tolist())
		# print 'next'
		ecc, a = fn.get_orbital_elements(particles, marr)



		# for i in [0,1]:
		# 	np.savetxt(save_file,np.c_[i, time, particles[i,0,0], particles[i,0,1], particles[i,0,2], \
		# 	particles[i,1,0], particles[i,1,1], particles[i,1,2], marr[i], v_kep], delimiter=',', fmt='%.10e')


		eccentricity.append(np.sqrt(sum(ecc)**2))
		semi_major_axis.append(a)
		vkep_list.append(v_kep)
		time_list.append(time)

		time += dt

		if a < 0.15:
			break

	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0


	return x_and_v, e_error, semi_major_axis, eccentricity, vkep_list, time_list


def hermite(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0
	iterations = 2
	x_and_v =[]
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


def plot_pos(particles, ax_range=1.2, tfinal=None, save=False, name='test'):
	'''
	Plots list (t, Np, 2, 3)
	'''
	xarr = np.array(particles)

	plt.plot(xarr[:, 0, 0, 0], xarr[:, 0, 0, 1], c='y')
	plt.plot(xarr[:, 1, 0, 0], xarr[:, 1, 0, 1], label='Test particle', c='indianred')

	# plt.scatter(1.,0, c="b", label='Initial postition')
	# plt.scatter(0,0, c="y", label='Sun')

	plt.legend()
	plt.xlabel('$x_{pos}$(au)')
	plt.ylabel('$y_{pos}$(au)')
	plt.xlim(-ax_range, ax_range)
	plt.ylim(-ax_range, ax_range)
	plt.title('Total integration time: %s years' %tfinal)
	plt.axis('equal')
	if save:
		plt.savefig('.png', transparant=True)
		plt.close()
	else:
		plt.show()


def plot_a(time, a, save=False):
	plt.plot(time, a, label='Test particle', c='indianred')
	plt.xlabel('Time (years)')
	plt.ylabel('Semi major axis (au)')
	plt.legend()
	if save:
		plt.savefig('5yrs_a.png', transparant=True)
		plt.close()
	else:
		plt.show()


def plot_error(timestep, tfinal):
	error1= []
	error2= []
	error3= []
	error4= []

	for dt in timestep:
		error1.append(euler   (dt,tfinal)  [1])
		error2.append(midpoint(dt,tfinal)  [1])
		error3.append(leapfrog(dt,tfinal,1)[1])
		error4.append(hermite (dt,tfinal)  [1])

	plt.plot(timestep, np.abs(error1), label = 'euler_forward')
	plt.plot(timestep, np.abs(error2), label = 'midpoint')
	plt.plot(timestep, np.abs(error3), label = 'leapfrog')
	plt.plot(timestep, np.abs(error4), label = 'hermite')

	plt.title('Error convergence')
	plt.legend()
	plt.xlabel("Stepsize [yr]")
	plt.ylabel("Error in energy")
	plt.xlim(1e-5,1e-2)
	plt.ylim(1e-13,1e-1)
	plt.yscale("Log")
	plt.xscale("Log")
	plt.axis('equal')
	plt.show()


def plot_vratio_vs_tau(tau, v_ratio,min=0, max=-1):

	plt.scatter(tau,np.average(v_ratio[int(min):int(max)]))

	tauvals = np.geomspace(1e-3, 1e3, num=100)
	plt.plot(tauvals, v_ratio_analytic(tauvals, 3))
	# gaussian(rad, r0=0.5, sigma=0.01, n=3)

	# plt.axvline(x=1, linestyle="dotted", c="indianred")
	plt.xscale("Log")
	plt.yscale("Log")
	plt.xlabel("$\\tau_{fric}$")
	plt.ylabel("$v_{r}/v_{K}$")
	plt.title("Radial drift for $n$=3")
	plt.legend()


def v_ratio_analytic(tau, n):
	"""
	Takes value for tau and returns ratio of vr and vk
	"""

	eta = n*0.05**2

	vrvk = -eta/(tau+(tau**(-1)))
	return np.abs(vrvk)



def load_files(dt, tfinal, tau):
	thing = np.zeros(2)
	particles = np.zeros([2,2,3])
	marr = np.zeros(2)
	v_kep = np.zeros(1)


	file_arr = np.loadtxt("object_time_x,v,m_dt=" + str(dt) + "_tfinal=" + str(tfinal) + "_tau=" +str(tau)+".csv",delimiter=',' )
	amount_of_steps = len(file_arr[:, 0])/2

	for line in file_arr:
		time = line[1]
		n = int(line[0])
		for dimension in range(3):
			particles[n, 0, dimension] = round(line[dimension + 2], 5)
			particles[n, 1, dimension] = round(line[dimension + 5], 5)
		marr[n] = line[-2]
		v_kep = line[-1]
		print marr
		print v_kep
		print particles, "\n\n"
	# return time, particles, v_kep

def vr_file(dt, tfinal, tau, save=True):
	pos_leapfrog, _, _,_, v_kep, time = leapfrog(dt, tfinal, tau, drag=True)

	dr = fn.calculate_vratio(dt, pos_leapfrog)

	if save:
		np.save("vratio_dt=" + str(dt)+"_tfinal=" + str(tfinal) + "_tau=" + str(tau),  [np.array(time[:-1]), dr/v_kep[:-1]] )

	return time[:-1], dr/v_kep[:-1]

def save_all_information(dt, tfinal, tau):
	'''save (all or some) information, savez saves the arrays seperatly, each accesible by keywords'''

	pos_leapfrog, _, a_leapfrog,_, v_kep, time = leapfrog(dt, tfinal, tau, drag=True)

	# np.save("pos_vkep_time_dt=" + str(dt) + "_tfinal=" + str(tfinal) + "_tau=" +str(tau), [np.array(pos_leapfrog), np.array(v_kep), np.array(time)])
	# np.savez("arrays_pos_vkep_time_dt=" + str(dt) + "_tfinal=" + str(tfinal) + "_tau=" +str(tau),  pos_leapfrog=pos_leapfrog, v_kep=v_kep, time=time )
	np.savez("eerste_pressure_bump_echt_met_bump_a_0.0148_dt=" + str(dt) + "_tfinal=" + str(tfinal) + "_tau=" +str(tau),  pos_leapfrog=pos_leapfrog, v_kep=v_kep, time=time, a_leapfrog=a_leapfrog)
	 # e_error=e_error, a_leapfrog=a_leapfrog, eccentricity=eccentricity)



def main():
	dt     = 0.0001
	tfinal = 40

	calc_step = 10
	omega_k = (2*np.pi)
	average_vratio = []

	tau_vals = np.geomspace(2e-3, 1e3, num=10)
	tau_vals[0] = 3.5e-3
 	tau_lijstje= list(tau_vals)
	tau_lijstje.append(1.)


	# tau_lijstje=[0.002,0.00859505945175,0.0369375234896, 0.158740105197,0.682190320772,2.93173318222,12.5992104989,54.1454816418,232.691816878,1000.0]

	# tau_we_need_high_tfinal = [3.5e-3, 0.008595059451754268,12.599210498948729,54.14548164181537,232.6918168776362,1000.0]
	tau_lijstje = [0.6821903207721965]


	for tau in tau_lijstje:
		print 'calculating', tau
		# save_all_information(dt, tfinal, tau)

		# leapfrog(dt, tfinal, tau, drag=True)

		# load_files(dt,tfinal,tau)
	 	# file_to_load = np.loadtxt( "object_time_x,v,m_dt=" + str(dt) + "_tfinal=" + str(tfinal) + "_tau=" +str(tau)+".csv",delimiter=',' )
		# print
		#		time, particles, v_kep = load_files(dt,tfinal,tau)


		# dr = fn.calculate_vratio(dt, particles)

		# plt.plot(time[:-1],  dr/v_kep[:-1] )


		# print probeer
		# print probeer1
		# print loaded_file
		# print np.loadtxt("hope_it_works_v_kep_dt=" + str(dt) + "_tfinal=" + str(tfinal) + "_tau=" +str(tau)+".dat")

		pvt_file = np.load("eerste_pressure_bump_echt_met_bump_a_0.0148_dt=" + str(dt) + "_tfinal=" + str(tfinal) + "_tau=" +str(tau)+".npz")


		print pvt_file.files

		dr = fn.calculate_vratio(dt, pvt_file['pos_leapfrog'])
		v_ratio = dr/pvt_file['v_kep'][:-1]

		# plt.plot( pvt_file['time'][:-1], v_ratio, label='{:0.2e}'.format(tau))
		# plot_vratio_vs_tau(tau, v_ratio, min=4/dt)


		# average_vratio.append(np.average(v_ratio[int(4/dt):]))
		# plt.plot(pvt_file['time'], pvt_file['a_leapfrog'], label='{:0.2e}'.format(tau))


	# print curve_fit(v_ratio_analytic, tau_lijstje, average_vratio)

		# '''use this to plot when files saved using np.save'''
		# pos_leapfrog, v_kep, time =\
		# np.load("pos_vkep_time_dt=" + str(dt) + "_tfinal=" + str(tfinal) + "_tau=" +str(tau)+".npy")
		# v_ratio = [ time[:-1], dr/v_kep[:-1] ]
		# plt.plot(v_ratio[0], v_ratio[1], label='{:0.2e}'.format(tau))



		plot_pos(pvt_file['pos_leapfrog'])
		# plt.show()

	plt.xlabel("time [yr]")
	plt.ylabel("vratio")
	# plt.xscale("Log")
	# plt.yscale("Log")
	# plt.ylabel("semi major axis $a$ [AU]")
	# # plt.title("vrad is (vji[0] - v_gas[0])* 0.1")
	plt.legend()
	plt.show()
	# # plt.savefig('semimajoraxis_extra_tau.png', transparant=True)
	# plt.close()

if __name__ == '__main__':
	main()
