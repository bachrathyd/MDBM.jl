# Unreleased

## ✨ New: vectorized (batched) evaluation

`MDBM_Problem(f, axes; vectorized = fv)` hands all new points of a stage (the
initial grid, each refinement, each neighbour check) to `fv(points)` in ONE call,
instead of evaluating `f` point by point in a threaded loop. Use it for GPU kernels
or any other batched evaluator. Fully backward compatible: without the keyword
nothing changes. See "Example 4" in the README and the `MDBM_Problem` docstring.

- `fv(points::Vector{<:Tuple})` returns the vector of function values (same type as `f`).
- `f` is still required (type detection, single-point calls) and must agree with `fv`.
- Requires memoization (the default); a constraint is still evaluated pointwise.
- Test: identical interpolated solution and evaluation cache to the scalar path.

Also documented: the default threaded batch evaluation needs a thread-safe `f`
(shared work arrays inside `f` cause data races that corrupt random points).

# MDBM.jl v0.2.2 Release Notes

This release introduces a major new feature: **error-based adaptive refinement**. This allows for more efficient and accurate discovery of solution manifolds by intelligently focusing computational effort on areas where the solution is changing most rapidly.

## ✨ New Features

### Error-Based Adaptive Refinement

The `solve!` function now supports error-based adaptive refinement, controlled by two new keyword arguments:

-   `refinementratio`: This parameter (between 0.0 and 1.0) controls the ratio of n-cubes to be refined based on their error. For example, a `refinementratio` of `0.3` will refine the 30% of n-cubes with the highest error.
-   `abstol`: This sets an absolute error tolerance for refinement. Only n-cubes with an error greater than this value will be refined.

These two parameters can be used together to achieve fine-grained control over the refinement process, leading to significant performance improvements and more accurate results, especially for complex problems.

#### Example

Here's how you can use the new error-based refinement feature:

```julia
using MDBM

# Define your problem...
# mymdbm = MDBM_Problem(...)

# Solve with error-based refinement:
# This will refine the 40% of n-cubes with the highest error in each iteration,
# but only if their error is above 1.0e-3.
solve!(mymdbm, 15, abstol=1.0e-3, refinementratio=0.4)
```

This adaptive approach avoids unnecessary refinement in smooth regions of the parameter space, concentrating the computational effort where it is most needed.

### High-Order Interpolation for Detailed Path Extraction

A new function `interpsubcubesolution!` has been introduced. It performs high-order interpolation within the solution n-cubes. This populates a tree structure (`posinterp`) in each n-cube, which can then be used with the `extract_paths` function to extract detailed geometric paths (e.g., curves in 3D problems).

### Dynamic Axis Extension

The `axesextend!` function has been improved to allow for more flexible extension of the parameter space. You can now provide a vector of new coordinates, and the function will automatically prepend or append them while maintaining the monotonicity of the axis.

## 🚀 Performance

The introduction of error-based adaptive refinement can lead to significant speedups, as fewer n-cubes are refined in each step, without sacrificing the accuracy of the solution.

We hope you enjoy these new features and improvements! As always, please feel free to open an issue on GitHub if you have any questions or feedback.
