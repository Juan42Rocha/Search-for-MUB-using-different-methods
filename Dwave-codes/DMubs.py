# 
# variacion de "Grafo, 2 entradas 2 salidas - Ising Model"
# de mau, aqui vamos a modificar para usar dwave,
#
#  Plan:
#   -comenzar bucando  matrices ortgonogales sin intentar optimizar
#   los las interacciones  
#

from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("agg")    # must select backend before imaporting pyplot
import matplotlib.pyplot as plt
from dimod import BinaryQuadraticModel
from dwave.system import LeapHybridSampler
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
import dwave.inspector
import neal
from time import time
import networkx as nx


def n_graph(n,fname,b,itr):
    """Returns a potential solution to the n-queens problem in a dictionary, where
    the keys are qubit id and value on or off

    Args:
        n: Number of queens to place.

        sampler: A binary quadratic model sampler. 
        Can use: 1) LeapHybridSampler
                 2) DWaveSampler - QPU
                 3) SimulatedAnnealingSampler
    """
    # read graph
    Jik = np.genfromtxt("../HamiltonianForIsingDwave/"+fname)#, delimiter=" ") 
    bqm = BinaryQuadraticModel({}, {}, 0, 'SPIN')
    bqm.offset = 0
    # iterations in the anealing
    #itr = 1#200    # mas interaciones mas tiempo?                    

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
    
    sampler = EmbeddingComposite(DWaveSampler())   
    sampleset = sampler.sample(bqm, num_reads=itr, label=f'{n} Ortogonalidad {start_time}')

    # Hybrid Solver
    # sampler = LeapHybridSampler()
    # sampleset = sampler.sample(bqm, label=f'{n}-q; {d}-d config')
    # print('quota_conversion_rate', sampler.properties['quota_conversion_rate'])

    # CPU
    #sampler = neal.SimulatedAnnealingSampler()
    #sampleset = sampler.sample(bqm, num_reads=itr)

    solver_time = f'finished in {time()-start_time} seconds'
    print(solver_time)
    # print('sampleset.info', sampleset.info)

    # ===============================================================================
    # ====== Save some infomation  ==================================================
    # ===============================================================================

    f = open(fname[0:-4]+"/Qs"+fname[0:-4]+"itr"+str(itr)+"state.txt", "w")
    minst = sampleset.record.sample
    for nqb in range(n):
        for i in range(itr):
            f.write(str(minst[i,nqb])+"   ")
        f.write("\n")
    f.close()

    #f = open("Css"+fname+"itr"+str(itr)+"b.txt", "w")
    #f.write(str(sampleset))
    #f.write(sampleset)
    #f.close()

    sample = sampleset.first.sample
    # print("Sample:\n", sample)
    f = open(fname[0:-4]+"/QoutGraph"+fname[0:-4]+"itr"+str(itr)+".txt", "a")
    f.write(f'{n} node graph\n')
    f.write(str(sample)+'\n')
    f.write(solver_time + '\n')
    f.write(str(sampleset.info)+'\n')
    f.write(str(sampleset.record.energy)+"\n")
    f.write(str(sampleset.record.num_occurrences)+"\n")
    f.close()

    # Inspect note: with a clasical solver the ispector doesnt work
    dwave.inspector.show(sampleset)
    return sampleset

"""
notas para el save  
    sampleset.record.sample   para array
                    .energy  para las energías 


"""

def is_valid_solution(n, solution):
    """Check energy

    Args:
        n: Number of queens in the problem.

        solution: A dictionary of qubits, key = qubit id, value = 0 or 1
    """
    sols = np.genfromtxt('sols.txt', delimiter=',')
    sol = np.array(list(solution.values()))
    print('Solution Found')
    print(sol.astype(float))
    print('Solution Hans')
    print(sols[0,:])
    print(sols[1,:])

    return np.array_equal(sol,sols[0,:]) or np.array_equal(sol,sols[1,:])

def plot_graph(n, sol):
    """Create a chessboard with queens using matplotlib. Image is saved
    in the root directory. Returns the image file name.
    """
    # Build graph for visualization
    # G = nx.Graph()
    G = nx.DiGraph()
    RelPad = np.genfromtxt('RelPad.txt', delimiter=',')
    # Jlr = np.genfromtxt('Jlr.txt', delimiter=',')
    # pos = {}

    dEdges = list()
    color_map = list()
    sol = list(sol.values())
    # print(sol)
    for i in range(n):
        # print(sol[i])
        G.add_node(i)
        if sol[i] == -1:
            color_map.append('r')
        else:
            color_map.append('b')
        # G.add_edge(i, RelPad[i,0]-1)
        # G.add_edge(i, RelPad[i,1]-1)
        dEdges.append((i,RelPad[i,0]-1))
        dEdges.append((i,RelPad[i,1]-1))
        # G.add_edges_from(i,(RelPad[i,0]-1))
        # G.add_edges_from(i,(RelPad[i,1]-1))
    G.add_edges_from(dEdges)


    # G = nx.Graph()

    # G.add_node(1)
    # G.add_node(2)
    # G.add_node(3)
    # G.add_node(4)
    # G.add_node(5)

    # G.add_edge(1, 2)
    # G.add_edge(1, 3)
    # G.add_edge(2, 3)
    # G.add_edge(2, 4)
    # G.add_edge(3, 4)
    # G.add_edge(3, 5)
    # G.add_edge(4, 5)
    # G.add_edge(4, 1)
    # G.add_edge(5, 1)
    # G.add_edge(5, 2)

    # pos = {
    #     0:(3,1),
    #     1:(4,3),
    #     2:(2,4),
    #     3:(0,3),
    #     4:(1,1),
    # }
    
    # nx.draw(G,pos=pos,with_labels=True,node_color="red",node_size=3000,
    #         font_color="white",font_size=20,font_family="DejaVu Sans",
    #         font_weight="bold",width=5)
    nx.draw(G, with_labels=True, node_color=color_map, pos=nx.shell_layout(G))
    # nx.draw(G, pos=nx.shell_layout(G))
    # plt.margins(0.2)
    # plt.show()

    file_name = 'grafo5.png'
    plt.savefig(file_name)

    return file_name

if __name__ == "__main__":

    d = 2
    N = 2
    k = 4
    Nqb = d*N*k**d

  # DHoHu2n1k2.dat   DHoHu2n1k4.dat      DHoHuHm2n2k4.dat
    # solo orth        ort en 4 fases      orth, mubs en 4 fases
    #    4                  26                  42         valores para b
    fname = "DHoHuHm2n2k4.dat"
    b = 42
    itr = 1
    print("Building Ising for Matrices {n} nodes.".format(n=Nqb))
    solutionset = n_graph(Nqb,fname,b, itr)

    minst = solutionset.record.sample
    print(minst)

 #   if is_valid_solution(n, solution):
 #       write = "Solution is valid."
 #   else:
 #        write = "Solution is invalid."
 #   
 #   print(write)
 #   f = open("outGraph.txt", "a")
 #   f.write(write + '\n\n')
 #   f.close()

 #   file_name = plot_graph(n, solution)
 #   print("Graph created. See: {}".format(file_name))
