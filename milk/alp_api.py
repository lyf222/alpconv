# coding: utf-8

import numpy as np
import ctypes
from ctypes import cdll
import matplotlib.pyplot as plt
import os

PATHOFTHISCODE = os.path.dirname(os.path.realpath(__file__))

def alp_api(gl, gb, smaxm, ma, ga, bfiledFlag, energies):
    """
    gl: galactic longitue, degree
    gb: galactic lattitude, degree
    smaxm: distance, kpc
    ma: neV
    ga: 10^11 GeV^-1
    energies: GeV
    """
    mylibc = cdll.LoadLibrary(PATHOFTHISCODE+"/alp/milk/alp.so")
    
    gl_c = ctypes.c_double(gl)
    gb_c = ctypes.c_double(gb)
    smaxm_c = ctypes.c_double(smaxm)
    ma_c = ctypes.c_double(ma)
    ga_c = ctypes.c_double(ga)
    bfiledFlag_c = ctypes.c_int(bfiledFlag)

    narray = len(energies)
    narray_c = ctypes.c_int(narray)
    array_c = (ctypes.c_double*narray)()
    energies_c = (ctypes.c_double*narray)()
    for i in range(narray):
        energies_c[i] = energies[i]
        array_c[i] = 0
    ptr_array_c = ctypes.pointer(array_c)
    ptr_energies_c = ctypes.pointer(energies_c)
    results = mylibc.alp(smaxm_c, gl_c, gb_c, ma_c, ga_c, bfiledFlag_c, ptr_energies_c, narray, ptr_array_c)
    prob = np.zeros_like(energies)
    for k in range(len(prob)):
        prob[k] = array_c[k]
    return prob

if __name__ == "__main__":
    #gl,gb,dist = 189.065,3.235,1.5
    #gl,gb,dist = 289.065,0.235,5.5
    #gl,gb,dist = 195.133,4.27009,0.25
    #gl,gb,dist = 75.2329, 0.112945,  5.0
    gl,gb,dist = 317.82 , -0.74, 8.5
    trueE = np.logspace(np.log10(1.e-6),np.log10(1e5),220)

    probs = alp_api(gl,gb,dist,20,1.5,1,trueE*1.e3)
    #plt.plot(trueE,probs,label='1')

    probs = alp_api(gl,gb,dist,1.,0.07,3,trueE)
    plt.plot(trueE,probs,label='10')

    probs = alp_api(gl,gb,dist,20,500,1,trueE)
    #plt.plot(trueE,probs,label='2')

    probs = alp_api(gl,gb,dist,20,2000,1,trueE)
    #plt.plot(trueE,probs,label='3')

    plt.legend()
    plt.xscale('log')
    plt.savefig('prob.png')
