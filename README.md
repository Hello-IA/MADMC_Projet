# Lorenz Non-Dominated Points for Multi-Objective Knapsack


## Requirements

- Python 3.7+
- NumPy
- Gurobi (with valid license)



## How to run:

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

