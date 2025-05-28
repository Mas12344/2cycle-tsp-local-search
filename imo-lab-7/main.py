import numpy as np
import random
import math
import matplotlib.pyplot as plt
from tqdm import tqdm
import ctypes
import os
import sys
import time

# g++ -shared -O3 -o local_search.dll local_search.cpp -static -static-libgcc -static-libstdc++ -std=c++20 -m64 -Wl,--subsystem,windows

class Instance(ctypes.Structure):
    _fields_ = [
        ("dist", (ctypes.c_int * 200) * 200)
    ]

class SolutionOut(ctypes.Structure):
    _fields_ = [
        ("cycle1", (ctypes.c_int * 100)),
        ("cycle2", (ctypes.c_int * 100)),
        ("value", (ctypes.c_int))
    ]

try:
    sa = ctypes.CDLL('./local_search.dll')
except Exception as e:
    print(f"DLL Loading Error: {e}")
    exit(1)


common_args = [
    Instance,
    ctypes.POINTER(SolutionOut),
]

common_args2 = [
    Instance,
    ctypes.POINTER(SolutionOut),
    ctypes.c_int
]

common_args3 = [
    Instance,
    ctypes.POINTER(SolutionOut),
    ctypes.c_int,
    ctypes.c_bool,
    ctypes.c_float,
]

hybrid_args = [
    Instance,
    ctypes.POINTER(SolutionOut),
    ctypes.c_int,
    ctypes.c_bool
]

mt_args = [
    Instance,
    ctypes.POINTER(SolutionOut),
    ctypes.c_int,
    ctypes.c_bool
]

sa.weighted_regret_heuristic.argtypes = common_args
sa.baseline_alg.argtypes = common_args
sa.multiple_start.argtypes = common_args2
sa.multiple_start.restype = ctypes.c_int
sa.iterated_alg.argtypes = common_args2
sa.iterated_alg.restype = ctypes.c_int
sa.destroy_repair_alg.argtypes = common_args3
sa.destroy_repair_alg.restype = ctypes.c_int
sa.hybrid_alg.argtypes = hybrid_args
sa.hybrid_alg.restype = ctypes.c_int
sa.mt_alg.argtypes = mt_args

def loadData(filename):
    file = open(filename, "r")
    lines = file.read().splitlines()

    size = 0
    points = []
    dataNow = False
    for line in lines:
        if(line == "EOF"): break
        if(dataNow):
            idx, x, y = line.split()
            points.append((float(x), float(y)))
        if(line == "NODE_COORD_SECTION"): dataNow = True
    return points

def euclidDistance(pointA, pointB):
    return round(math.dist(pointA, pointB), 0)

def genInst(points):
    inst = Instance()
    for i in range(len(points)):
        for j in range(len(points)):
            if(i == j): inst.dist[i][j] = 0
            else:
                dist = euclidDistance(points[i], points[j])
                inst.dist[i][j] = int(dist)

    return inst

def genDistMatrix(points):
    distMatrix = []
    for i in range(len(points)):
        distMatrix.append([])
        for j in range(len(points)):
            if(i == j): distMatrix[i].append(0)
            else:
                dist = euclidDistance(points[i], points[j])
                distMatrix[i].append(round(dist))
    return distMatrix


def random_init():
    cycle = list(range(200))
    random.shuffle(cycle)

    cycle1 = cycle[:100]
    cycle2 = cycle[100:]

    sol = SolutionOut()
    sol.value = -1
    for i,c in enumerate(cycle1):
        sol.cycle1[i] = c
    for i,c in enumerate(cycle2):
        sol.cycle2[i] = c
    return sol


def print_results_value(results, count):
    results_list = results.items()
    kroA = []
    kroB = []
    for k,v in results_list:
        best_val = min([s.value for s,t in v])
        worst_val = max([s.value for s,t in v])

        mean_val = sum([s.value for s,t in v]) / count
        r = [k[8:], f"{mean_val} ({best_val} - {worst_val})"]
        if k.startswith("kroA"):
            kroA.append(r)
        else:
            kroB.append(r)
    
    print("kroA200")
    for w in kroA:
        print(f"\t{w[0]}: {w[1]}")
    print("kroB200")
    for w in kroB:
        print(f"\t{w[0]}: {w[1]}")

def save_best(results, kroa, krob):
    def save_plot(solution, points, name):
        cycle1, cycle2 = [*solution.cycle1], [*solution.cycle2]
        plt.clf()
        for cycle, color in zip([cycle1, cycle2], ['r', 'b']):
            x = [points[i][0] for i in cycle + [cycle[0]]]
            y = [points[i][1] for i in cycle + [cycle[0]]]
            plt.plot(x, y, marker='o', color=color)
        plt.title(name)
        plt.axis('off')
        plt.savefig(f'{name}.png')

    for k,v in results.items():
        best_val = float("inf")
        best_sol = None
        for s,t in v:
            if s.value < best_val:
                best_val = s.value
                best_sol = s

        if k.startswith("kroA"):
            save_plot(best_sol, kroa, k)
        else:
            save_plot(best_sol, krob, k)


def print_results_time(results):
    results_list = results.items()
    kroA = []
    kroB = []
    for k,v in results_list:
        best_val = min([t for s,t in v])
        worst_val = max([t for s,t in v])
        mean_val = sum([t for s,t in v])/len(v)
        r = [k[8:], f"{round(mean_val, 2)} ({round(best_val, 2)} - {round(worst_val, 2)})"]
        if k.startswith("kroA"):
            kroA.append(r)
        else:
            kroB.append(r)
    
    print("kroA200")
    for w in kroA:
        print(f"\t{w[0]}: {w[1]}")
    print("kroB200")
    for w in kroB:
        print(f"\t{w[0]}: {w[1]}")


def main():
    results = {}
    kroA, kroB = loadData("kroA200.tsp"), loadData("kroB200.tsp")
    k = 10

    for instance_name, data in zip(["kroA200", "kroB200"], [kroA, kroB]):
        inst = genInst(data)
        distMatrix = genDistMatrix(data)

        for algorithm in ["mt_alg"]:
            id_name = f"{instance_name}_{algorithm}"
            print(id_name)
            tmp_res = []

            multiple_times = []
            timed = 0
            for i in tqdm(range(k)):
                output = random_init()

                args = inst, ctypes.byref(output)
                
                start = time.perf_counter()
                if algorithm == "weighted_regret_heuristic":
                    sa.weighted_regret_heuristic(*args)
                elif algorithm == "hybrid_alg":
                    iters = sa.hybrid_alg(*args, 6800, True)
                elif algorithm == "hybrid_alg (no local search)":
                    iters = sa.hybrid_alg(*args, 6800, False)
                elif algorithm == "mt_alg":
                    iters = sa.mt_alg(*args, 6800, True)
                elif algorithm == "mt_alg (no local search)":
                    iters = sa.mt_alg(*args, 6800, False)
                elif algorithm == "baseline":
                    sa.baseline_alg(*args)
                elif algorithm == "cache":
                    sa.cache_alg(*args)
                elif algorithm == "multiple_start":
                    timed = sa.multiple_start(*args, 200)
                elif algorithm == "iterated_alg":
                    iters = sa.iterated_alg(*args, 6800)
                elif algorithm == "destroy_repair_alg (no local search) 0.3":
                    sa.destroy_repair_alg(*args, 6800, False, 0.3)
                elif algorithm == "destroy_repair_alg 0.3":
                    sa.destroy_repair_alg(*args, 6800, True, 0.3)
                elif algorithm == "destroy_repair_alg (no local search) 0.5":
                    sa.destroy_repair_alg(*args, 6800, False, 0.5)
                elif algorithm == "destroy_repair_alg 0.5":
                    sa.destroy_repair_alg(*args, 6800, True, 0.5)
                elif algorithm == "destroy_repair_alg (no local search)":
                    sa.destroy_repair_alg(*args, 6800, False, 0.6)
                elif algorithm == "destroy_repair_alg":
                    sa.destroy_repair_alg(*args, 6800, True, 0.6)
                
                elapsed = time.perf_counter() - start

                tmp_res.append([output, elapsed*1000])
                multiple_times.append(timed/1000)
            results.update({id_name: tmp_res})

    print(f"\n{multiple_times=}")
    print(f"avg = {sum(multiple_times) / len(multiple_times)}")
    print("="*15)
    print_results_value(results, k)
    print("time [ms]")
    print_results_time(results)
    save_best(results, kroA, kroB)

if __name__ == '__main__':
    main()