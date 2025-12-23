from read_file import readFile
import numpy as np
w = np.zeros(16,dtype=int)
v = np.zeros((16,2),dtype=int)
filename = "2KP200-TA-0_test.dat"
capacity, weights, values = readFile(filename,w,v)

def dominates(a, b):
    return all(x >= y for x, y in zip(a.objectives, b.objectives)) and any(x > y for x, y in zip(a.objectives, b.objectives))

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
    print(v)
    t = [[[] for _ in range(capacity + 1)] for _ in range(v + 1)]

    t[0][0] = [Solution((0,)*v)]  
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

def Lorenz_otp(solution):
    opt = []
    for s in solution:
        obj = s.objectives[:]
        obj = tuple(sorted(obj))
        sum = 0
        lorenz = []
        for o in obj:
            sum += o
            lorenz.append(sum)
        opt.append(
                        Solution(
                            objectives=tuple(lorenz),
                            prev=s.prev,
                            taken=s.taken
                        )
                    )
    return opt

DP = pareto_dinamique(capacity, weights, values)

lorenz_solution = Lorenz_otp(DP[len(values)][capacity])
for sol in lorenz_solution:
    s = sol
    choix = []
    while s.prev is not None:
        choix.append(s.taken)
        s = s.prev
    choix.reverse()
    print("vecteur solution", sol.prev.objectives, end="")
    print("Lorenz vecteur solution ", sol.objectives, end="")
    print(choix)
