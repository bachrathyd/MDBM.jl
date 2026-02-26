# MDBM.jl: Multi-Dimensional Bisection Method

**MDBM.jl** is a robust and efficient Julia package for finding roots of implicit non-linear equation systems. It is specifically designed to determine high-dimensional submanifolds (points, curves, surfaces, etc.) of solutions.

---

## 🚀 Quick Start

```julia
using MDBM

# 1. Define your function (roots of f(x,y) = 0)
function f(x, y)
    return x^2 + y^2 - 4.0
end

# 2. Define the parameter space (Axes)
ax1 = Axis(-5.0:1.0:5.0, "x")
ax2 = Axis(-5.0:1.0:5.0, "y")

# 3. Create and solve the problem
mdbm = MDBM_Problem(f, [ax1, ax2])
solve!(mdbm, 4) # 4 iterations of refinement

# 4. Get the results
x_sol, y_sol = getinterpolatedsolution(mdbm)
```

---

## 🛠 Core API

### `Axis`
Defines a dimension in the parameter space.
- `Axis(range, name, is_periodic=false)`
- Support for `StepRangeLen`, `Vector`, and even `Vector{Vector{T}}` (Vector-valued axes).

### `MDBM_Problem`
The central container for the problem.
- `MDBM_Problem(f, axes; constraint=c, memoization=true, contour_level_fc=0.0)`
- `f`: Function returning a scalar or a Tuple/Vector of values.
- `constraint`: A function that must return `> 0` for the solution to be valid.

### `solve!`
Iteratively refines the solution.
- `solve!(mdbm, iterations; interpolationorder=1, checkneighbourNum=1)`

---

## 🌟 Advanced Features

### 1. Periodic Boundaries
Useful for angular coordinates ($[0, 2\pi]$).
```julia
ax_theta = Axis(0:0.1:2pi, "theta", true) # true marks it as periodic
```

### 2. Sub-cube Interpolation & Path Extraction
Go beyond the grid resolution to extract smooth paths or surfaces.
```julia
interpsubcubesolution!(mdbm)
paths = extract_paths(mdbm) # Returns individual "branches" of the solution
```

### 3. Vector-Valued Axes
Search along a path in a higher-dimensional space as if it were a single dimension.
```julia
# A 1D axis where each point is actually a 2D coordinate
ax_v = Axis([[0,0], [1,1], [2,0]], "path") 
```

### 4. Dynamic Axis Extension
Expand your search space while keeping previous results.
```julia
axesextend!(mdbm, 1, append=[6.0, 7.0, 8.0])
```

### 5. Custom Contour Levels
Find solutions for $f(x) = c$ instead of $f(x) = 0$.
```julia
mdbm = MDBM_Problem(f, axes, contour_level_fc=1.5)
```

---

## 📚 Examples Gallery

| File | Description | Key Feature |
| :--- | :--- | :--- |
| `base_par4_codim2.jl` | 4 Parameters, 2 Equations | High-dimension, Path extraction |
| `periodical_coordinates.jl` | Mapping a sphere | Periodic axes |
| `vector_axis.jl` | Searching on triangle edges | Vector-valued ticks |
| `eigenvalue_based_problem.jl` | Stability analysis | Non-linear root finding |
| `milling_stability_chart.jl` | Mechanical engineering | Real-world application |
| `interactive_neighbour_check.jl` | Interactive Visualization | Error-based refinement |

---

## 💡 Performance Tips

1. **Initial Resolution**: Use a coarse initial grid. MDBM will refine only where the solution exists.
2. **Memoization**: Always keep memoization on for expensive functions.
3. **Multi-threading**: Set `ENV["JULIA_NUM_THREADS"]` to speed up the `solve!` process.
4. **Neighbor Check**: If your solution manifold is disconnected or very thin, use `checkneighbourNum > 1`.

---

## 📄 Citation

If you use MDBM.jl in your research, please cite:
> Bachrathy, Dániel. "Multi-Dimensional Bisection Method." *Periodica Polytechnica Mechanical Engineering* (2018). [Link](https://pp.bme.hu/me/article/view/1236/640)
