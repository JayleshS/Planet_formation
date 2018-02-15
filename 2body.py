import numpy as np
import pars
import functions as fn
import matplotlib.pyplot as plt

# xarr, varr, marr = fn.init_2body(0)
# xarr, varr, marr = fn.init_2body(0)
# etot0 = fn.e_tot(xarr, varr, marr)

tfinal = 0.002*pars.yr
dt = 0.001*pars.yr
# etot0 = fn.e_tot(xarr, varr, marr)


x1, y1, x2, y2 = ([] for i in range(4))

def euler(dt, tfinal):
	xarr, varr, marr = fn.init_2body(0)
	etot0 = fn.e_tot(xarr, varr, marr)
	time = 0
	while time < tfinal:
		acc = fn.forces(xarr, varr, marr)
		xarr += varr*dt
		varr += acc*dt
		x1.append(xarr[0, 0])
		y1.append(xarr[1, 0])
		x2.append(xarr[0, 1])
		y2.append(xarr[1, 1])
		time += dt
	etot1 = fn.e_tot(xarr, varr, marr)
	e_error = (etot1 - etot0) / etot0
	return x1, y1, x2, y2, e_error

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

		x1.append(xarr[0, 0])
		y1.append(xarr[1, 0])
		x2.append(xarr[0, 1])
		y2.append(xarr[1, 1])
		# coords.append(xarr)

		time += dt

	etot1 = fn.e_tot(xarr, varr, marr)
	e_error = (etot1 - etot0) / etot0
	print x1
	print y1
	print x2
	print y2 
	return x1, y1, x2, y2, e_error, np.array(coords)
	# print dingems
	return x1, y1, x2, y2, e_error, np.array(dingems)

def leapfrog(dt, tfinal):
	xarr, varr, marr = fn.init_2body(0)
	etot0 = fn.e_tot(xarr, varr, marr)
	dingems = [xarr]
	time = 0
	# print dingems
	while time < tfinal:
		acc = fn.forces(xarr, varr, marr)


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

	# plot(euler(xarr, varr, marr)[:-1])
	# for time in [0.1*pars.yr, 0.01*pars.yr, 0.001*pars.yr]:
	x1_mid, y1_mid, x2_mid, y2_mid, e_error_mid, ding = midpoint(dt, tfinal)
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
