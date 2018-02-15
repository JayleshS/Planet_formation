import numpy as np
import pars
import matplotlib.pyplot as plt


def init_2body(ecc):
    xarr = np.zeros((3, pars.Np))
    varr = np.zeros((3, pars.Np))
    marr = np.zeros(pars.Np)
    marr[0] = pars.mSun
    marr[1] = pars.mEarth

    # vKep = np.sqrt((pars.mEarth+pars.mSun)/pars.a1)
    # vKep = sqrt((pars.muEarth+pars.muSun)/sqrt((2*pars.au**2)))
    vKep = np.sqrt(pars.gN*(pars.mSun + pars.mEarth) / pars.au)

    xarr[:, 1] = [pars.au * (1+ecc), 0., 0.]
    # xarr[:, 1] = [pars.au * (1+ecc), 0., .5*pars.au*(1+ecc)]

    varr[:, 1] = [0., vKep * np.sqrt((1-ecc)/(1+ecc)), 0.]

    return xarr, varr, marr
'''
def forces(xarr, varr, marr):
    acc = np.zeros((3,pars.Np))
    for i in range(pars.Np):
        for j in range(i+1, pars.Np):
            rji = xarr[:,j] - xarr[:,i]
            r2 = sum(rji**2)
            r1 = np.sqrt(r2)
            r3 = r1*r2

            force = pars.gN*rji/r3
            acc[:,i] += force*marr[j]
            acc[:,j] -= force*marr[i]

    return acc



def e_tot(xarr, varr, marr):
    Ekin = 0
    Epot = 0
    for i in range(pars.Np):
        vi2 = sum(varr[:, i]**2)
        Ekin += 0.5*marr[i]*vi2
        for j in range(i+1, pars.Np):
            rji = xarr[:, j] - xarr[:, i]
            r1 = np.sqrt(sum(rji**2))
            Epot += -(pars.gN*marr[i]*marr[j])/r1
    Etot = Ekin + Epot
    return Etot
'''
def forces(xarr, varr, marr):
    acc = np.zeros((3,pars.Np))
    for i in range(pars.Np):
        for j in range(i+1, pars.Np):
            rji = xarr[:,j] - xarr[:,i]

            acc[:,i] += pars.gN*rji/sum(rji**2)**(3/2) * marr[j]
            acc[:,j] -= pars.gN*rji/sum(rji**2)**(3/2) * marr[i]

    return acc


def e_tot(xarr, varr, marr):
    Ekin = 0
    Epot = 0
    for i in range(pars.Np):
        Ekin += 0.5*marr[i]*sum(varr[:, i]**2)
        for j in range(i+1, pars.Np):
            Epot += -(pars.gN*marr[i]*marr[j]) / (np.sqrt(sum((xarr[:, j] - xarr[:, i]))**2))
    Etot = Ekin + Epot
    return Etot
