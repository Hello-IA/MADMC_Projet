# Multi-Objective Knapsack Problem - Lorenz Efficiency Solver

This project implements algorithms for finding Lorenz-efficient solutions to the multi-objective knapsack problem using Gurobi optimization.

## Overview

The project contains implementations of:
- **P1**: Finds a single Lorenz-efficient solution using Ordered Weighted Averaging (OWA)
- **PL**: Iteratively finds multiple Lorenz-nondominated solutions with dominance cuts
- **Lorenz Vector Enumeration**: Enumerates all solutions with the same Lorenz vector (for p=2 objectives)

## Requirements

### Dependencies

- Python 3.7+
- NumPy
- Gurobi Optimizer (with valid license)

### Installation

1. Install required packages:
```bash
pip install numpy gurobipy
```

2. Ensure you have a valid Gurobi license. For academic use, you can obtain a free license from [Gurobi's website](https://www.gurobi.com/academia/academic-program-and-licenses/).

## File Structure

- `P1Model.py` - Main implementation file containing:
  - `P1()` - Single Lorenz-efficient solution finder
  - `PL()` - Iterative Lorenz solution finder with dominance cuts
  - `iterative_PL()` - Wrapper for iterative execution
  - `enumerate_same_lorenz_p2()` - Enumerate all x solutions for a given Lorenz vector (p=2)
  - `add_no_good_cut()` - Helper function to add no-good cuts

- `read_file.py` - File reader for benchmark instances (required)

## Input Data Format

The solver expects data files in the `.dat` format with:
- Knapsack capacity
- Item weights
- Item values (multi-objective)

Example: `2KP200-TA-0_test.dat`

## Usage

### Basic Usage - Running P1Model.py

To run the main script with default settings:

```bash
python P1Model.py
```

This will:
1. Read the test instance `2KP200-TA-0_test.dat`
2. Solve for multiple omega weight vectors: [1.01, 1], [1.5, 1], [2, 1], [5, 1], [10, 1]
3. Find all Lorenz-efficient solutions for each omega
4. Enumerate different x solutions for each Lorenz vector
5. Print results showing OWA objective values and Lorenz vectors

### Using Individual Functions

#### 1. Finding a Single Lorenz-Efficient Solution (P1)

```python
import numpy as np
from P1Model import P1
from read_file import readFile

# Load data
w = np.zeros(16, dtype=int)
v = np.zeros((16, 2), dtype=int)
capacity, weights, values = readFile("2KP200-TA-0_test.dat", w, v)

# Define omega weights (must be strictly decreasing and positive)
omega = [2.0, 1.0]

# Solve P1
solution = P1(weights, values, capacity, omega, verbose=True)

print("Optimal x:", solution["x"])
print("Objective values y:", solution["y"])
print("Lorenz vector L:", solution["L"])
print("OWA objective:", solution["obj"])
```

#### 2. Finding Multiple Lorenz-Efficient Solutions (Iterative PL)

```python
from P1Model import iterative_PL

omega = [2.0, 1.0]

# Find all Lorenz-efficient solutions
solutions, Lorenz_vectors = iterative_PL(
    weights, values, capacity, omega,
    verbose=True,
    max_iters=None  # No limit on iterations
)

print(f"Found {len(solutions)} Lorenz-efficient solutions")
for i, sol in enumerate(solutions):
    print(f"Solution {i}: L={sol['L']}, OWA={sol['obj']}")
```

#### 3. Enumerating Solutions for a Fixed Lorenz Vector (p=2)

```python
from P1Model import enumerate_same_lorenz_p2

# Use a Lorenz vector from previous solutions
L_star = solutions[0]["L"]

# Find all x with same Lorenz vector
same_L_solutions = enumerate_same_lorenz_p2(
    weights, values, capacity, L_star,
    verbose=True,
    max_solutions=100  # Limit number of solutions
)

print(f"Found {len(same_L_solutions)} different x solutions with L={L_star}")
```

## Parameters

### P1 / PL Parameters

- `weights`: Array of item weights (1D numpy array)
- `values`: Array of item values for each objective (2D numpy array, shape n×p)
- `capacity`: Knapsack capacity (float)
- `omega`: Weight vector for OWA (strictly decreasing, positive values)
- `time_limit`: Optional time limit in seconds (default: None)
- `mip_gap`: Optional MIP optimality gap (default: None)
- `verbose`: Print solver output (default: True)

### iterative_PL Additional Parameters

- `max_iters`: Maximum number of iterations (default: None - no limit)

### enumerate_same_lorenz_p2 Additional Parameters

- `max_solutions`: Maximum number of solutions to enumerate (default: None - no limit)

## Understanding the Output

Each solution dictionary contains:
- `"x"`: Binary decision vector (which items to include)
- `"y"`: Objective values for each criterion
- `"L"`: Lorenz vector (cumulative sorted objective values)
- `"obj"`: OWA objective value
- `"status"`: Gurobi solver status

### Lorenz Vector Interpretation

For p=2 objectives with values y₁ and y₂:
- L[0] = min(y₁, y₂) - worst objective value
- L[1] = y₁ + y₂ - sum of both objectives

A solution is Lorenz-dominated if there exists another solution with a componentwise better L vector.

## Example Output

```
omega: [2, 1]
Found 5 Lorenz vectors (in Lorenz space).
  [0] OWA=15234.000000  L=[7500. 15234.]
  [1] OWA=15100.000000  L=[7600. 15200.]
  [2] OWA=14980.000000  L=[7650. 15150.]
  ...
  
  Lorenz #0: L*=[7500. 15234.] -> 3 different x solutions
     y= [7500. 7734.]
     y= [7734. 7500.]
     y= [7617. 7617.]
```

## Troubleshooting

### Common Issues

1. **"No module named 'gurobipy'"**
   - Install Gurobi: `pip install gurobipy`
   - Verify license is properly configured

2. **"No module named 'read_file'"**
   - Ensure `read_file.py` is in the same directory as `P1Model.py`

3. **Infeasible model**
   - Check that capacity is sufficient
   - Verify omega values are strictly decreasing and positive
   - Ensure input data is correctly formatted

4. **Slow performance**
   - Set `time_limit` parameter
   - Set `mip_gap` to allow suboptimal solutions (e.g., 0.01 for 1% gap)
   - Reduce `max_iters` or `max_solutions` for enumeration

## Algorithm Details

### P1 - Initial Lorenz-Efficient Solution
Uses a linearized OWA formulation to maximize Lorenz dominance. Any optimal solution is guaranteed to be Lorenz-nondominated.

### PL - Iterative Lorenz Solution Finder
Adds dominance cuts to exclude previously found Lorenz vectors and find new ones. The algorithm terminates when no improving solution exists.

### Mathematical Formulation
- Uses dual variables r_k and b_{k,i} to linearize the Lorenz dominance relation
- Objective: maximize Σ_k λ_k [(k+1)r_k - Σ_i b_{k,i}]
- Where λ is derived from omega: λ_k = ω_k - ω_{k+1} (and λ_p = ω_p)

## References

This implementation is based on optimization techniques for multi-objective knapsack problems using Lorenz efficiency and OWA operators.

## License

Academic use only.
