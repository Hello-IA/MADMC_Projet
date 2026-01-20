import pdb
from turtle import pd
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from read_file import readFile
import numpy as np

def P1(weights, values, capacity, omega, time_limit=None, mip_gap=None, verbose=True):
    
    """Find one Lorenz-efficient solution of the multi-objective knapsack problem. 
       Any optimal solution of P1 is Lorenz-nondominated
    """
    weights = np.asarray(weights, dtype=float)
    values  = np.asarray(values, dtype=float)

    n = len(weights)
    p = values.shape[1]
    
    assert values.shape[0] == n, "values must be shape (n,p)"
    assert len(omega) == p, "omega must have length p"
    assert all(omega[i] > omega[i+1] for i in range(p-1)) and all(w > 0 for w in omega), \
        "omega must be strictly decreasing and positive"

    # lambda from omega: (ω1-ω2, ω2-ω3, ..., ω_{p-1}-ω_p, ω_p)
    lam = np.zeros(p, dtype=float)
    for k in range(p-1):
        lam[k] = omega[k] - omega[k+1]
    lam[p-1] = omega[p-1]
    assert np.all(lam > 0), "lambda must be strictly positive"

    m = gp.Model("P1_Lorenz_OWA")

    if not verbose:
        m.Params.OutputFlag = 0
    if time_limit is not None:
        m.Params.TimeLimit = time_limit
    if mip_gap is not None:
        m.Params.MIPGap = mip_gap

    # x_j ∈ {0,1}
    x = m.addVars(n, vtype=GRB.BINARY, name="x")

    y = m.addVars(p, lb=0.0, vtype=GRB.CONTINUOUS, name="y")

    # r_k ∈ R and b_{k,i} ≥ 0
    r = m.addVars(p, lb=0.0, vtype=GRB.CONTINUOUS, name="r")      # lb=0 ok if values >= 0
    b = m.addVars(p, p, lb=0.0, vtype=GRB.CONTINUOUS, name="b")   # b[k,i] = b_i^k
    
   


    # Capacity: sum w_j x_j <= W
    m.addConstr(gp.quicksum(weights[j] * x[j] for j in range(n)) <= capacity, name="capacity")

    #  y_i = sum_j v_{j,i} x_j 
    for i in range(p):
        m.addConstr(
            y[i] == gp.quicksum(values[j, i] * x[j] for j in range(n)),
            name=f"y_def[{i}]"
        )


    # r_k - b_{k,i} <= y_i   for all k,i
    for k in range(p):
        for i in range(p):
            m.addConstr(
                r[k] - b[k, i] <= y[i],
                name=f"lorenz_dual[{k},{i}]"
            )
            
    # max sum_k lam[k] * ( (k+1)*r[k] - sum_i b[k,i] )
    obj = gp.LinExpr()
    for k in range(p):
        obj += lam[k] * ((k + 1) * r[k] - gp.quicksum(b[k, i] for i in range(p)))

    m.setObjective(obj, GRB.MAXIMIZE)


    m.optimize()

    # if m.Status not in [GRB.OPTIMAL, GRB.TIME_LIMIT]:
    #     # For debugging: INFEASIBLE would indicate a modeling error here (P1 should be feasible)
    #     return None

    x_sol = np.array([int(round(x[j].X)) for j in range(n)], dtype=int)
    y_sol = np.array([y[i].X for i in range(p)], dtype=float)


    y_sorted = np.sort(y_sol)
    L_sol = np.cumsum(y_sorted)

    return {
        "x": x_sol,
        "y": y_sol,
        "L": L_sol,
        "obj": m.ObjVal if m.SolCount > 0 else None,
        "status": m.Status,
    }

    
    


def PL(weights, values, capacity, omega, prev_Ls, time_limit=None, mip_gap=None, verbose=True):
    """now forbid solutions whose Lorenz vector is dominated by a previously found one.
    """
    weights = np.asarray(weights, dtype=float)
    values  = np.asarray(values, dtype=float)

    n = len(weights)
    p = values.shape[1]
    l = len(prev_Ls)

    # lambda from omega: (ω1-ω2, ω2-ω3, ..., ω_{p-1}-ω_p, ω_p)
    lam = np.zeros(p, dtype=float)
    for k in range(p-1):
        lam[k] = omega[k] - omega[k+1]
    lam[p-1] = omega[p-1]
    assert np.all(lam > 0), "lambda must be strictly positive"

    m = gp.Model("PL_Lorenz_OWA")

    if not verbose:
        m.Params.OutputFlag = 0
    if time_limit is not None:
        m.Params.TimeLimit = time_limit
    if mip_gap is not None:
        m.Params.MIPGap = mip_gap

    # x_j ∈ {0,1}
    x = m.addVars(n, vtype=GRB.BINARY, name="x")

    y = m.addVars(p, lb=0.0, vtype=GRB.CONTINUOUS, name="y")

    # r_k ∈ R and b_{k,i} ≥ 0
    r = m.addVars(p, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="r")
     # lb=0 ok if values >= 0
    b = m.addVars(p, p, lb=0.0, vtype=GRB.CONTINUOUS, name="b")   # b[k,i] = b_i^k
    
    z = m.addVars(l, p, vtype=GRB.BINARY, name="z")


    # Capacity: sum w_j x_j <= W
    m.addConstr(gp.quicksum(weights[j] * x[j] for j in range(n)) <= capacity, name="capacity")

    #  y_i = sum_j v_{j,i} x_j 
    for i in range(p):
        m.addConstr(
            y[i] == gp.quicksum(values[j, i] * x[j] for j in range(n)),
            name=f"y_def[{i}]"
        )

    # Lorenz dual constraints: r_k - b_{k,i} <= y_i
    for k in range(p):
        for i in range(p):
            m.addConstr(r[k] - b[k, i] <= y[i],
                        name=f"lorenz_dual[{k},{i}]")

    
    # Lorenz-component expressions: T_k = (k+1) r_k - sum_i b_{k,i}
    T = [ (k+1)*r[k] - gp.quicksum(b[k,i] for i in range(p)) for k in range(p) ]

    # Dominance cuts vs each previously found Lorenz vector prev_Ls[s]
    # For each s: choose at least one k to improve, and enforce improvement if chosen.
    for s in range(l):
        m.addConstr(gp.quicksum(z[s, k] for k in range(p)) >= 1,
                    name=f"one_improvement[{s}]")
        for k in range(p):
            # If z[s,k] = 1 => T_k >= prev_Ls[s][k] + 1
            # If z[s,k] = 0 => constraint is T_k >= 0 (since RHS=0)
            m.addConstr(T[k] >= (prev_Ls[s][k] + 1) * z[s, k],
                        name=f"improve[{s},{k}]")

            
    # max sum_k lam[k] * ( (k+1)*r[k] - sum_i b[k,i] )
    obj = gp.LinExpr()
    for k in range(p):
        obj += lam[k] * ((k + 1) * r[k] - gp.quicksum(b[k, i] for i in range(p)))

    m.setObjective(obj, GRB.MAXIMIZE)

    m.optimize()
    
    x_sol = np.array([int(round(x[j].X)) for j in range(n)], dtype=int)
    y_sol = np.array([y[i].X for i in range(p)], dtype=float)


    y_sorted = np.sort(y_sol)
    L_sol = np.cumsum(y_sorted)

    return {
        "x": x_sol,
        "y": y_sol,
        "L": L_sol,
        "obj": m.ObjVal if m.SolCount > 0 else None,
        "status": m.Status,
    }




w = np.zeros(16,dtype=int)
v = np.zeros((16,2),dtype=int)
filename = "2KP200-TA-0_test.dat"
capacity, weights, values = readFile(filename,w,v)
objectives = values.shape[1]
omegas = [[a, 1] for a in [1.01, 1.5, 2, 5, 10]]

for omega in omegas:
    sol1 = P1(weights, values, capacity, omega, verbose=False)
    prev_Ls = [sol1["L"]]
    sol2 = PL(weights, values, capacity, omega, prev_Ls, verbose=False)

    print("omega:", omega, 
          "P1 OWA:", sol1["obj"],
          "P1 L:", sol1["L"],
          "PL OWA:", sol2["obj"] if sol2["obj"] else None)
"""we see that in the PL's solution, the worst-case is slightly worse but the total is better than in P1's solution."""