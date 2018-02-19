import numpy as np
import pars
import functions_matrix_form as fn
import matplotlib.pyplot as plt

# xarr, varr, marr = fn.init_2body(0)
# xarr, varr, marr = fn.init_2body(0)
# etot0 = fn.e_tot(xarr, varr, marr)

tfinal = 0.002*pars.yr
dt = 0.001*pars.yr
# etot0 = fn.e_tot(xarr, varr, marr)

x_and_v =[]

def euler(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0
	while time < tfinal:
		acc = fn.forces(particles, marr)
		particles[:,0,:] += particles[:,1,:]*dt
		particles[:,1,:] += acc*dt

		x_and_v.append(particles.tolist())

		time += dt
	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0
	return x_and_v, e_error

euler(1.,10)



def midpoint(dt, tfinal):
	particles, marr = fn.init_2body(0)
	etot0 = fn.e_tot(particles, marr)
	time = 0

	while time < tfinal:
		acc = fn.forces(particles, marr)
		particles_mid = np.zeros((pars.Np, 2, 3))
		particles_mid[:, 0, :] = particles[:,0,:] + particles[:,1,:] * dt / 2
		particles_mid[:, 1, :] = particles[:,1,:] + acc * dt / 2

		amid = fn.forces(particles_mid, marr)
		particles[:, 0, :] += particles_mid[:, 1, :] * dt
		particles[:,1,:] += amid * dt

		x_and_v.append(particles.tolist())

		time += dt

	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	print 'dag'

	return particles, e_error

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
		partciles[:,1,:] += acc* dt/2

		time += dt
	etot1 = fn.e_tot(particles, marr)
	e_error = (etot1 - etot0) / etot0

	return particles

def plot(x1_val, y1_val, x2_val, y2_val):
	# plt.plot(x1_val, y1_val)
	# plt.plot(x2_val, y2_val)
	# print x2_val
	# plt.show()
	pass

def plottest(xarr):
	# plt.plot(xarr[:,0], xarr[:,1])
	# plt.plot(xarr[:, 0, 0], xarr[:, 0, 1])
	# plt.plot(xarr[:, 1, 0], xarr[:, 1, 1])
	# plt.show()
	# print xarr[:, 1, 0]
	# print xarr
	pass







def main():


	# x1_euler, y1_euler, x2_euler, y2_euler, e_error_euler = euler(xarr, varr, marr)
	# plot(x1_euler, y1_euler, x2_euler, y2_euler)
	# print e_error_euler

	# plot(euler(positions, marr)[:-1])
	# for time in [0.1*pars.yr, 0.01*pars.yr, 0.001*pars.yr]:
	'''x_mid,v_mid,  e_error_mid = midpoint(dt, tfinal)'''
	# plot(x1_mid, y1_mid, x2_mid, y2_mid)
	# plottest(ding)

	# print ding
	# plot(x1_mid, y1_mid, x2_mid, y2_mid)
	# print e_error_euler
		# print e_error_mid
	# print ding[:,1]
	# gr = [3]
	# f = 0
	# while f < 5:
	# 	print gr
	# 	gr.append(f*2)
	# 	f += 1
	# print gr

if __name__ == '__main__':
	main()
