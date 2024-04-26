from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("agg")    # must select backend before importing pyplot
import matplotlib.pyplot as plt
from dimod import BinaryQuadraticModel
from dwave.system import LeapHybridSampler
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
import dwave.inspector
import neal
from time import time
import networkx as nx


d = 2
N = 1
K = 2
Nqb = d*N*K**d

# DHoHu2n1k2.dat   DHoHu2n1k4.dat      DHoHuHm2n2k4.dat
# solo orth        ort en 4 fases      orth, mubs en 4 fases
#    4                  26                  42         valores para b
fname = "DHoHu2n1k2.dat"
b = 4
itr = 3

Jik = np.genfromtxt("../HamiltonianForIsingDwave/"+fname)#, delimiter=" ") 
bqm = BinaryQuadraticModel({}, {}, 0, 'SPIN')
bqm.offset = 0

n = Nqb
# ===============================================================================
# ======= Build a graph with from file fname ====================================
# ===============================================================================
fl = 0
for i in range(n):
    bqm.add_variable(i, b)
    for k in range(n):
        if Jik[i,k] > 1e-14:
            bqm.add_interaction(i, k, Jik[i,k])    
            fl = 1

print(bqm)

if fl == 0:
    print("Error en Jik")
    exit()

# ===============================================================================
# ====== Selecto solver method ==================================================
# ===============================================================================

print('solver started...')
start_time = time()
    
# QPU select the sampler
# sampler = EmbeddingComposite(DWaveSampler(solver=dict(topology__type='pegasus')))
# sampler = EmbeddingComposite(DWaveSampler(solver=dict(topology__type='zephyr')))
    
#sampler = EmbeddingComposite(DWaveSampler())   
#sampleset = sampler.sample(bqm, num_reads=itr, label=f'{n} Ortogonalidad {start_time}')

# Hybrid Solver
# sampler = LeapHybridSampler()
# sampleset = sampler.sample(bqm, label=f'{n}-q; {d}-d config')
# print('quota_conversion_rate', sampler.properties['quota_conversion_rate'])

# CPU
sampler = neal.SimulatedAnnealingSampler()
sampleset = sampler.sample(bqm, num_reads=itr)

solver_time = f'finished in {time()-start_time} seconds'
print(solver_time)
# print('sampleset.info', sampleset.info)


