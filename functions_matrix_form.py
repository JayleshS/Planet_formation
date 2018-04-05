import numpy as np
import pars
import matplotlib.pyplot as plt


def init_2body(ecc):
    particles = np.zeros((pars.Np,2,3))
    marr = np.zeros(pars.Np)
    marr[0] = pars.mSun
    marr[1] = pars.mEarth

    v_kep = np.sqrt(pars.gN*(pars.mSun + pars.mEarth) / pars.au)

    particles[1,0,:] = [pars.au * (1+ecc), 0., 0.]  # Initial position of 2nd particle
    particles[1,1,:] = [0., v_kep * np.sqrt((1-ecc)/(1+ecc)), 0.]  # Initial velocity of 2nd particle

    return particles, marr

def forces(particles, marr):
    acc = np.zeros((pars.Np,3))
    for i in range(pars.Np):
        for j in range(i+1, pars.Np):
            rji = particles[j, 0, :] - particles[i, 0, :]
            acc[i, :] += pars.gN * rji / (np.sqrt(sum(rji**2)) * sum(rji**2)) * marr[j]
            acc[j, :] -= pars.gN * rji / (np.sqrt(sum(rji**2)) * sum(rji**2)) * marr[i]

    return acc


def forces_hermite(particles, marr):
    acc = np.zeros((pars.Np, 3))
    jer = np.zeros((pars.Np, 3))
    for i in range(pars.Np):
        for j in range(i+1, pars.Np):
            rji = particles[j, 0, :] - particles[i, 0, :]
            vji = particles[j, 1, :] - particles[i, 1, :]

            r2 = sum(rji**2)
            r1 = np.sqrt(r2)
            r3 = r1 * r2
            rv = sum(rji * vji)
            rv /= r2

            acc[i, :] += pars.gN * (rji / r3) * marr[j]
            acc[j, :] -= pars.gN * (rji / r3) * marr[i]

            jerk = (vji - 3 * rv * rji) / r3
            jer[i, :] += pars.gN * jerk * marr[j]
            jer[i, :] -= pars.gN * jerk * marr[i]

    return acc, jer

def forces_migration(particles, marr, tau):
    acc = np.zeros((pars.Np,3))
    rji = particles[1, 0, :] - particles[0, 0, :]
    vji = particles[1, 1, :] - particles[0, 1, :]

    # Convert to cylindrical coordinates
    rad = np.sqrt(rji[0]**2 + rji[1]**2)
    theta = np.arctan2(rji[1], rji[0])

    v_kep  = np.sqrt(pars.gN*(marr[0] + marr[1]) / rad)
    v_head = 0.004 * v_kep

    '''9975, 9987, 0.004'''

    t_stop = tau / (v_kep / rad)

    v_gas    = np.zeros(3)
    v_gas[0] = -np.sin(theta) * (v_kep - v_head) #+ np.cos(theta) * (vji[0] - v_gas[0]) * 0.1
    v_gas[1] =  np.cos(theta) * (v_kep - v_head) #+ np.sin(theta) * (vji[1] - v_gas[1]) * 0.1

    acc[1, :] = - (vji - v_gas) / t_stop
    return acc, v_kep



def forces_total(particles, marr,tau):
    acc_tot = np.zeros((pars.Np, 3))

    acc_grav = forces(particles, marr)
    acc_mig, v_kep = forces_migration(particles, marr, tau)
    acc_tot = acc_grav + acc_mig

    return acc_tot, v_kep



def e_tot(particles, marr):
    Ekin = 0
    Epot = 0
    for i in range(pars.Np):
        Ekin += 0.5*marr[i]*sum(particles[i,1,:]**2)
        for j in range(i+1, pars.Np):
            Epot += -(pars.gN*marr[i]*marr[j]) / (np.sqrt(sum((particles[j,0,:] - particles[i,0,:])**2)))
    Etot = Ekin + Epot
    return Etot


def get_orbital_elements(particles, marr):
    ang = np.zeros((pars.Np, 3))
    n = np.zeros(3)

    rji = particles[1, 0, :] - particles[0, 0, :]
    vji = particles[1, 1, :] - particles[0, 1, :]
    r1 = np.sqrt(sum(rji**2))

    ang[0, :] += np.cross(rji, vji)
    ang[1, :] += np.cross(rji, vji)
    lz = ang[1, 2]
    inc = np.arccos(lz / np.sqrt(sum(ang[1, :]**2)))

    e = (np.cross(vji, ang[1, :])) / (pars.gN*sum(marr)) - (rji/r1)

    a = (sum(ang[1, :]**2)) / (pars.gN*sum(marr)*(1-sum(e**2)))

    return e, a
