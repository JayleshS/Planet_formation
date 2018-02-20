import numpy as np
import pars
import functions as fn
import matplotlib.pyplot as plt

# xarr, varr, marr = fn.init_2body(0)
# xarr, varr, marr = fn.init_2body(0)
# etot0 = fn.e_tot(xarr, varr, marr)

tfinal = 1.1*pars.yr
dt = 0.001*pars.yr
# etot0 = fn.e_tot(xarr, varr, marr)

# x = []
# v = []



def euler(dt, tfinal):
	xarr, varr, marr = fn.init_2body(0)
	etot0 = fn.e_tot(xarr, varr, marr)
	time = 0
	x1 = []
	x2 = []
	y1 = []
	y2 = []
	while time < tfinal:
		acc = fn.forces(xarr, varr, marr)
		xarr += varr*dt
		varr += acc*dt

		x1.append(xarr[0, 0])
		y1.append(xarr[1, 0])
		x2.append(xarr[0, 1])
		y2.append(xarr[1, 1])

		# x.append(list(xarr))
		# v.append(list(varr))



		time += dt
	etot1 = fn.e_tot(xarr, varr, marr)
	e_error = (etot1 - etot0) / etot0
	return x1, y1, x2, y2, e_error

# print 'hi'
# euler(5.,10)


def midpoint(dt, tfinal):
	xarr, varr, marr = fn.init_2body(0)
	etot0 = fn.e_tot(xarr, varr, marr)
	coords = [xarr]
	time = 0

	while time < tfinal:
		acc = fn.forces(xarr, varr, marr)
		xmid = xarr + varr * dt / 2
		vmid = varr + acc * dt / 2
		amid = fn.forces(xmid, vmid, marr)
		xarr += vmid * dt
		varr += amid * dt

		x.append(list(xarr))
		v.append(list(varr))
		# coords.append(xarr)

		time += dt

	etot1 = fn.e_tot(xarr, varr, marr)
	e_error = (etot1 - etot0) / etot0

	# print dingems
	return x, v, e_error

def leapfrog(dt, tfinal):
	xarr, varr, marr = fn.init_2body(0)
	etot0 = fn.e_tot(xarr, varr, marr)
	time = 0
	# print dingems
	while time < tfinal:
		acc   = fn.forces(xarr, varr, marr)
		varr += acc* dt/2
		xarr += varr*dt
		acc   = fn.forces(xarr, varr, marr)
		varr += acc* dt/2

		time += dt
	etot1 = fn.e_tot(xarr, varr, marr)
	e_error = (etot1 - etot0) / etot0

	return xarr, varr

def plot(x1_val, y1_val, x2_val, y2_val):
	plt.plot(x1_val, y1_val)
	plt.plot(x2_val, y2_val)
	# print x2_val
	plt.show()
	# pass

def plottest(xarr):
	# plt.plot(xarr[:,0], xarr[:,1])
	# plt.plot(xarr[:, 0, 0], xarr[:, 0, 1])
	# plt.plot(xarr[:, 1, 0], xarr[:, 1, 1])
	# plt.show()
	# print xarr[:, 1, 0]
	# print xarr
	pass



def main():


	x1_euler, y1_euler, x2_euler, y2_euler, e_error_euler = euler(dt, tfinal)
	plot(x1_euler, y1_euler, x2_euler, y2_euler)
	# print e_error_euler

	# plot(euler(xarr, varr, marr)[:-1])

if __name__ == '__main__':
	main()
