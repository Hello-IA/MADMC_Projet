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
    # r = m.addVars(p, lb=0.0, vtype=GRB.CONTINUOUS, name="r")
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

    m.Params.DualReductions = 0
    
    # U = float(np.max(np.sum(values, axis=0)))  # safe upper bound on any y_i
    # for k in range(p):
    #     m.addConstr(r[k] <= U, name=f"r_ub[{k}]")
    m.optimize()


    if m.Status == GRB.INFEASIBLE:
        print("infeasible")
        return {"status": GRB.INFEASIBLE, "obj": None, "x": None, "y": None, "L": None}
    if m.Status == GRB.UNBOUNDED:
        print("unbounded")
        raise RuntimeError("PL is unbounded (modeling issue).")
    if m.Status == GRB.TIME_LIMIT and m.SolCount == 0:
        print("time limit reached with no solution")
        raise RuntimeError("TIME_LIMIT reached but no feasible solution found.")
    if m.SolCount == 0:
        raise RuntimeError(f"No solution available. Status={m.Status}")

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

def iterative_PL(weights, values, capacity, omega, time_limit=None, mip_gap=None, verbose=False, max_iters=None):

    sols = []

    # Phase 1: P1
    sol = P1(weights, values, capacity, omega,
             time_limit=time_limit, mip_gap=mip_gap, verbose=verbose)

    if sol["status"] == GRB.INFEASIBLE or sol["obj"] is None:
        raise RuntimeError(f"P1 failed (status={sol['status']}).")

    sols.append(sol)
    Ls = [sol["L"]]

    it = 0
    while True:
        if max_iters is not None and it >= max_iters:
            break

        sol_pl = PL(weights, values, capacity, omega, Ls,
                    time_limit=time_limit, mip_gap=mip_gap, verbose=verbose)


        if sol_pl["status"] == GRB.INFEASIBLE:
            break

        sols.append(sol_pl)
        Ls.append(sol_pl["L"])
        it += 1

        if verbose:
            print(f"[it={it}] status={sol_pl['status']} obj={sol_pl['obj']} L={sol_pl['L']}")

    return sols, Ls


def add_no_good_cut(m, x, x_sol, name="nogood"):
    expr = gp.LinExpr()
    for j, val in enumerate(x_sol):
        if int(val) == 1:
            expr += (1 - x[j])
        else:
            expr += x[j]
    m.addConstr(expr >= 1, name=f"{name}_{m.NumConstrs}")


def enumerate_same_lorenz_p2(weights, values, capacity, L_star,
                             time_limit=None, mip_gap=None, verbose=False, max_solutions=None):
    """
    Enumerate all different x such that Lorenz(y)=L*
    Lorenz(y) = (min(y1,y2), y1+y2).
    """
    weights = np.asarray(weights, dtype=float)
    values  = np.asarray(values, dtype=float)

    n = len(weights)
    p = values.shape[1]
    assert p == 2, "This helper is for p=2 only."

    L1, L2 = float(L_star[0]), float(L_star[1])

    m = gp.Model("FIX_L_ENUM_X_P2")
    if not verbose:
        m.Params.OutputFlag = 0
    if time_limit is not None:
        m.Params.TimeLimit = float(time_limit)
    if mip_gap is not None:
        m.Params.MIPGap = float(mip_gap)

    x = m.addVars(n, vtype=GRB.BINARY, name="x")
    y = m.addVars(2, lb=0.0, vtype=GRB.CONTINUOUS, name="y")

    # capacity
    m.addConstr(gp.quicksum(weights[j] * x[j] for j in range(n)) <= capacity, name="capacity")

    # y definition
    for i in range(2):
        m.addConstr(y[i] == gp.quicksum(values[j, i] * x[j] for j in range(n)), name=f"y_def[{i}]")

    # Fix L2 = y1 + y2
    m.addConstr(y[0] + y[1] == L2, name="fix_sum")

    # Enforce min(y0,y1) = L1:
    # y0 >= L1 and y1 >= L1
    m.addConstr(y[0] >= L1, name="min_lb0")
    m.addConstr(y[1] >= L1, name="min_lb1")

    # And force at least one equals L1 (via <= with big-M and a binary)
    # If t=1 -> y0 <= L1 ; if t=0 -> y1 <= L1
    t = m.addVar(vtype=GRB.BINARY, name="t_min")

    # Safe upper bound for y components:
    U = float(np.max(np.sum(values, axis=0)))  # very safe
    M = U  # big-M

    m.addConstr(y[0] <= L1 + M * (1 - t), name="min_eq0")
    m.addConstr(y[1] <= L1 + M * (t),     name="min_eq1")

    # Objective irrelevant (feasibility enumeration). Use 0.
    m.setObjective(0.0, GRB.MAXIMIZE)

    sols = []
    seen = set()

    while True:
        if max_solutions is not None and len(sols) >= max_solutions:
            break

        m.optimize()

        if m.Status == GRB.INFEASIBLE:
            break
        if m.SolCount == 0:
            raise RuntimeError(f"No feasible solution for fixed Lorenz {L_star}. status={m.Status}")

        x_sol = np.array([int(round(x[j].X)) for j in range(n)], dtype=int)
        xk = tuple(x_sol.tolist())
        if xk in seen:
            add_no_good_cut(m, x, x_sol, name="nogood_dup")
            continue
        seen.add(xk)

        y_sol = np.array([y[i].X for i in range(2)], dtype=float)

        sols.append({"x": x_sol, "y": y_sol, "L": np.array([L1, L2])})

        add_no_good_cut(m, x, x_sol, name="nogood")

    return sols




if __name__ == "__main__":
    w = np.zeros(16, dtype=int)
    v = np.zeros((16, 2), dtype=int)

    filename = "2KP200-TA-0_test.dat"
    capacity, weights, values = readFile(filename, w, v)

    objectives = values.shape[1]
    omegas = [[a, 1] for a in [1.01, 1.5, 2, 5, 10]]

    for omega in omegas:
        sols, Ls = iterative_PL(
            weights, values, capacity, omega,
            verbose=False
        )
        for idx, L_star in enumerate(Ls):
            sameL = enumerate_same_lorenz_p2(weights, values, capacity, L_star, verbose=False)
            print(f"  Lorenz #{idx}: L*={L_star} -> {len(sameL)} different x solutions")
            for s in sameL[:5]:
                print("     y=", s["y"])

        print("\n==============================")
        print("omega:", omega)
        print("Found", len(sols), "Lorenz vectors (in Lorenz space).")

        for t, sol in enumerate(sols):
            print(f"  [{t}] OWA={sol['obj']:.6f}  L={sol['L']}")


        sol1 = sols[0]
        sol_last = sols[-1]
        print("P1 OWA:", sol1["obj"], "P1 L:", sol1["L"])
        print("Last OWA:", sol_last["obj"], "Last L:", sol_last["L"])




# w = np.zeros(16,dtype=int)
# v = np.zeros((16,2),dtype=int)
# filename = "2KP200-TA-0_test.dat"
# capacity, weights, values = readFile(filename,w,v)
# objectives = values.shape[1]
# omegas = [[a, 1] for a in [1.01, 1.5, 2, 5, 10]]

# for omega in omegas:
#     sol1 = P1(weights, values, capacity, omega, verbose=False)
#     prev_Ls = [sol1["L"]]
#     sol2 = PL(weights, values, capacity, omega, prev_Ls, verbose=False)

#     print("omega:", omega, 
#           "P1 OWA:", sol1["obj"],
#           "P1 L:", sol1["L"],
#           "PL OWA:", sol2["obj"] if sol2["obj"] else None)
# """we see that in the PL's solution, the worst-case is slightly worse but the total is better than in P1's solution."""