# Lorenz Efficiency Solver for Multi-Objective Knapsack

Find Lorenz-efficient solutions to multi-objective knapsack problems using Gurobi.

## Requirements

- Python 3.7+
- NumPy
- Gurobi (with valid license)

```bash
pip install numpy gurobipy
```

Get a free academic Gurobi license at: https://www.gurobi.com/academia/

## Quick Start

Run the main script:

```bash
python P1Model.py
```

This will solve the test instance and display all Lorenz-efficient solutions.

## Main Functions

### 1. Dynamic Programming

Finds one Lorenz-efficient solution.

```
python lorenz_non_domines_indirecte.py

```

### 2. iterative_PL and enumerate_same_lorenz_p2 - Alternative Solutions

Finds all Lorenz-efficient solutions.

```
python P1Model.py
```


## Output

Each solution contains:
- `x` - Items to select (binary array)
- `y` - Objective values
- `L` - Lorenz vector (sorted cumulative values)
- `obj` - OWA objective value

## Parameters

**Required:**
- `weights` - Item weights (numpy array)
- `values` - Item values per objective (2D array)
- `capacity` - Knapsack capacity
- `omega` - Weight vector (strictly decreasing)

