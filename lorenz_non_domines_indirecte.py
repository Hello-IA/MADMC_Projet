from read_file import readFile
import numpy as np
w = np.zeros(16,dtype=int)
v = np.zeros((16,2),dtype=int)
filename = "2KP200-TA-0_test.dat"
capacity, weights, values = readFile(filename,w,v)

def owa_score(obj, omega):
    y_sorted = sorted(obj)  # ascending
    return sum(w * v for w, v in zip(omega, y_sorted))

def omega_default(p):
    return [p - i for i in range(p)]  # [p, p-1, ..., 1]

def dominates(a, b):
    return all(x >= y for x, y in zip(a.objectives, b.objectives)) and any(x > y for x, y in zip(a.objectives, b.objectives))

def dominates_(a, b):
    return all(x >= y for x, y in zip(a, b)) and any(x > y for x, y in zip(a, b))

def pareto_filter(solutions):
    pareto = []
    for s in solutions:
        dominated = False
        to_remove = []
        for p in pareto:
            if dominates(p, s):
                dominated = True
                break
            if dominates(s, p):
                to_remove.append(p)
        if not dominated:
            for r in to_remove:
                pareto.remove(r)
            pareto.append(s)
    return pareto


class Solution:
    def __init__(self, objectives, prev=None, taken=False):
        self.objectives = tuple(objectives)
        self.prev = prev
        self.taken = taken



def pareto_dinamique(capacity, weights, values):
    v = len(values)
    p = values.shape[1]  # number of objectives
    # print(v)
    t = [[[] for _ in range(capacity + 1)] for _ in range(v + 1)]

    t[0][0] = [Solution((0,) * p)]  #dimention of objectives use 
    
    for i in range(1, v+1): 
        wi = weights[i - 1]
        vi = values[i - 1]  

        for w in range(capacity + 1):

            solutions = []

            for s in t[i - 1][w]:
                solutions.append(
                    Solution(
                        objectives=s.objectives,
                        prev=s,
                        taken=False
                    )
                )


            if wi <= w:
                for s in t[i - 1][w - wi]:
                    new_obj = tuple(a + b for a, b in zip(s.objectives, vi))
                    solutions.append(
                        Solution(
                            objectives=new_obj,
                            prev=s,
                            taken=True
                        )
                    )

            t[i][w] = pareto_filter(solutions)
    return t

def Lorenz_vec(solution):
    srt = sorted(solution)   
    acc = 0
    opt = []
    for s in srt:
        acc += s
        opt.append(acc)
    return tuple(opt)

def lorenz_transform(pareto_solutions):
    return [Solution(Lorenz_vec(s.objectives), prev=s, taken=False) for s in pareto_solutions]





print("---------Part1---------")
# DP[len(values)][capacity]: this is pareto front for full capacity, not =< capacity
n = len(weights)
DP = pareto_dinamique(capacity, weights, values)
all_solutions = [s for w in range(capacity + 1) for s in DP[n][w]]
pareto_front = pareto_filter(all_solutions)
lorenz_solutions = lorenz_transform(pareto_front)

objectives = values.shape[1]
omega = omega_default(objectives)

# ... after computing lorenz_solutions ...

for sol in lorenz_solutions:
    base = sol.prev
    choix = []
    s = base
    while s.prev is not None:
        choix.append(s.taken)
        s = s.prev
    choix.reverse()

    owa = owa_score(base.objectives, omega)

    print(
        "Objective:", base.objectives,
        "Lorenz:", sol.objectives,
        "OWA:", owa,   
        "Choices:", choix
    )

print("---------Part1:counterexample---------")
# Part 1 counterexample:
solution1 = lorenz_solutions[2]
solution2 = lorenz_solutions[4]
print("solution 1 Objective:", solution1.prev.objectives, "Lorenz:", solution1.objectives)
print("solution 2 Objective:", solution2.prev.objectives, "Lorenz:", solution2.objectives)

if dominates_(solution1.objectives, solution2.objectives):
    print("Solution 1 Lorenz-dominates Solution 2")
else:
    print("Solution 1 does not Lorenz-dominate Solution 2")
    
added_item = (0, 876)
solution1_added = tuple(
    a + b for a, b in zip(solution1.prev.objectives, added_item))
solution2_added = tuple(
    a + b for a, b in zip(solution2.prev.objectives, added_item)
)
print("After adding (0,876) to both:")
print("new Solution 1 objectives:", solution1_added, " Lorenz:", Lorenz_vec(solution1_added))
print("new Solution 2 objectives:", solution2_added, " Lorenz:", Lorenz_vec(solution2_added))

    
if dominates_(solution1_added, solution2_added):
    print("Solution 1 Lorenz-dominates Solution 2")
else:
    print("Solution 1 does not Lorenz-dominate Solution 2")