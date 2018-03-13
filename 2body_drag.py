import numpy as np
import pars
import functions_matrix_form as fn
import matplotlib.pyplot as plt
# import matplotlib as mpl

plt.rcParams['axes.linewidth'] = 2   # Create thicker lines around plots
plt.rcParams['lines.linewidth'] = 2  # Create thicker lines in plots

# xarr, varr, marr = fn.init_2body(0)
# xarr, varr, marr = fn.init_2body(0)
# etot0 = fn.e_tot(xarr, varr, marr)


# etot0 = fn.e_tot(xarr, varr, marr)

x_and_v =[]
eccentricity_save = []
semi_major_axis_save = []


def leapfrog(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0

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
		acc = grav + drag

		particles[1,1,:] += acc* dt/2
		particles[1,0,:] += particles[1,1,:]*dt

		grav   = fn.forces(particles, marr)[1]
		drag = fn.forces_migration(particles, marr)
		acc = grav + drag*1e1

		particles[1,1,:] += acc* dt/2

		eccentricity, semi_major_axis = fn.get_orbital_elements(particles, marr)


		eccentricity_save.   append(eccentricity   .tolist())
		semi_major_axis_save.append(semi_major_axis.tolist())
		x_and_v.append(particles.tolist())

		time += dt

	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return x_and_v, e_error, eccentricity_save, semi_major_axis_save





def plot(particles):
	'''
	Plots list (t, Np, 2, 3)
	'''
	xarr = np.array(particles)
	for planet in range(pars.Np):
		plt.plot(xarr[:, planet, 0, 0], xarr[:, planet, 0, 1])
	plt.plot(xarr[:, planet, 0, 0][0:1100], xarr[:, planet, 0, 1][0:1100], c='r', label="Initial orbit")

	plt.legend()
	plt.xlabel('$x_{pos}$[cm]')
	plt.ylabel('$y_{pos}$[cm]')
	plt.show()

def main():

 	# plot( leapfrog_drag( 0.01*pars.yr, 120.0 * pars.yr)[0] )


	plt.plot(leapfrog_drag( 0.001*pars.yr, 1.0 * pars.yr)[3], label="a")
	plt.legend()
	plt.show()


if __name__ == '__main__':
	main()
