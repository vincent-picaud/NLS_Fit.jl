# Do not `using NLS_Solver`

The reason is that, as `NLS_Fit` already embeds the `NLS_Solver.jl`
package and you must use this embedded package version. This is the
recommended approach to avoid version conflicts...

Contritely, you just have to add a `NLS_Solver` prefix when you need
some of the `NLS_Solver.jl` package functionalities. By example:
```@example
conf = NLS_Solver.LevenbergMarquardt_Conf()
```

See [`NLS_ForwardDiff_From_Model2Fit`](@ref) for furthers details.



