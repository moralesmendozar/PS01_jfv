# -----------------------------------------------------------------------------
# COMPUTATIONAL ECONOMICS
# Homework 1, Fall 2019
# Prof. Jesus Fernandez Villaverde
# Due date: December 2nd, 2019
#
# Author: Juan Castellanos Silván
# Date: November 22nd, 2019
# Julia, Version 1.2.0
# ------------------------------------------------------------------------------

# ---------
# Packages
# ---------
using Parameters, Plots

# --------
# Modules
# -------
include("steady_state.jl")
using .exercise2

include("fixed_grid.jl")
using .exercise3

#include("egm_generalized.jl")
#using .exercise4

include("accelerator.jl")
using .exercise6

include("multi_grid.jl")
using .exercise7

# ---------------
# 1. Calibration
# ---------------
villaverde_economy = @with_kw (
                        # a) Preferences
                        β = 0.96, # time discount factor
                        θ = 0.5,  # elasticity of substitution between goods

                        # b) Technology
                        α = 0.33, # capital share
                        δ = 0.1, # depreciation

                        # c) Stochastic processes
                        # production func. 1
                        vGridZ = [-0.0673, -0.0336, 0, 0.0336, 0.0673],
                        mTransitionZ = [0.9727 0.0273 0 0 0;
                                       0.0041 0.9806 0.0153 0 0;
                                       0 0.0082 0.9836 0.0082 0;
                                       0 0 0.0153 0.9806 0.0041;
                                       0 0 0 0.0273 0.9727],
                        # production func. 2
                        vGridA = [0.9, 1, 1.1],
                        mTransitionA = [0.9 0.1 0;
                                        0.05 0.9 0.05;
                                        0 0.1 0.9]
                        )

villaverdeEconomy = villaverde_economy()

# ------------------------------
# 2. Deterministic steady state
# ------------------------------
@unpack α, δ, θ = villaverdeEconomy

# Solving system of equations
initialGuess = [0.5, 0.5, 1]
results = exercise2.ss_solve(villaverdeEconomy, initialGuess)
laborOneSteadyState = results.zero[1]
laborTwoSteadyState = results.zero[2]
capitalSteadyState = results.zero[3]

# From the resource constraints
consumptionOneSteadyState = capitalSteadyState^(α) * laborOneSteadyState^(1-α) - δ * capitalSteadyState
consumptionTwoSteadyState = laborTwoSteadyState


villaverde_SteadyState = @with_kw (
                        capitalSteadyState = capitalSteadyState,

                        laborOneSteadyState = laborOneSteadyState,

                        consumptionOneSteadyState = consumptionOneSteadyState,

                        laborTwoSteadyState = laborTwoSteadyState,

                        consumptionTwoSteadyState = consumptionTwoSteadyState,

                        utilitySteadyState = consumptionOneSteadyState^(θ) * consumptionTwoSteadyState^(1-θ) - 0.5*(laborOneSteadyState + laborTwoSteadyState)^2
)

villaverdeSteadyState = villaverde_SteadyState()

println(" Deterministic steady state ")
println(" -------------------------- ")
println(" K_ss = ", round(capitalSteadyState, digits= 5), "  c_ss^1 = ",round(consumptionOneSteadyState, digits=5), "  l_ss^1 = ", round(laborOneSteadyState, digits=5), "  c_ss^2 = ", round(consumptionTwoSteadyState, digits=5), "  l_ss^2 = ", round(laborTwoSteadyState, digits=5))
println(" ")

# --------------
# 3. Fixed grid
# --------------
nGridCapital = 250

@time tValueFunctionFixedGrid, tCapitalPolicyFixedGrid, vGridCapital, vMaxDifferenceFixedGrid = exercise3.fixed_grid(villaverdeEconomy, villaverdeSteadyState, nGridCapital)

# Graphs
# 3.1.
pValueFunctionFixedGrid = plot(vGridCapital, tValueFunctionFixedGrid[:,1,1], label = "z_1, A_1", xlabel = "capital", legend=:bottomright)
plot!(vGridCapital, tValueFunctionFixedGrid[:,end,1], label = "z_5, A_1")
plot!(vGridCapital, tValueFunctionFixedGrid[:,1,end], label = "z_1, A_3")
plot!(vGridCapital, tValueFunctionFixedGrid[:,end,end], label = "z_5, A_3")

savefig(pValueFunctionFixedGrid, "/Users/Castesil/Documents/EUI/Year II - PENN/Computational Economics/Homework 1/LaTeX/pValueFunctionFixedGrid.pdf")

# 3.2.
pCapitalPolicyFixedGrid = plot(vGridCapital, tCapitalPolicyFixedGrid[:,1,1], label = "z_1, A_1", xlabel = "capital", legend=:bottomright)
plot!(vGridCapital, tCapitalPolicyFixedGrid[:,end,1],label = "z_5, A_1")
plot!(vGridCapital, tCapitalPolicyFixedGrid[:,1,end], label = "z_1, A_5")
plot!(vGridCapital, tCapitalPolicyFixedGrid[:,end,end], label = "z_5, A_3")
plot!(vGridCapital, vGridCapital, color=:black, linestyle=:dash, label= "45 degrees")

savefig(pCapitalPolicyFixedGrid, "/Users/Castesil/Documents/EUI/Year II - PENN/Computational Economics/Homework 1/LaTeX/pCapitalPolicyFixedGrid.pdf")

# 3.3.
pConvergenceFixedGrid = plot(vMaxDifferenceFixedGrid[2:end], ylabel = "error", xlabel = "number of iterations", legend=:none)
savefig(pConvergenceFixedGrid, "/Users/Castesil/Documents/EUI/Year II - PENN/Computational Economics/Homework 1/LaTeX/pConvergenceFixedGrid.pdf")

# --------------------------
# 4. Endogenous Grid Method
# --------------------------
nGridCapital = 60

# 4.1. hold labor constant at ss
tValueTilde = exercise4steadystate.egm_labor_ss(villaverdeEconomy, villaverdeSteadyState, nGridCapital)

# 4.2. recover policies


# ------------------------
# 5. Comparison of grids
# ------------------------

# See pdf file with solutions

# ----------------
# 6. Accelerator
# ----------------
nGridCapital = 250

@time tValueFunctionAccelerator, tCapitalPolicyAccelerator, vGridCapital, vMaxDifferenceAccelerator = exercise6.accelerator(villaverdeEconomy, villaverdeSteadyState, nGridCapital)

# Graphs
# 6.1.
pCapitalPolicyAccelerator = plot(vGridCapital, tCapitalPolicyAccelerator[:,1,1], label = "z_1, A_1", xlabel = "capital", legend=:bottomright)
plot!(vGridCapital, tCapitalPolicyAccelerator[:,end,1],label = "z_5, A_1")
plot!(vGridCapital, tCapitalPolicyAccelerator[:,1,end], label = "z_1, A_5")
plot!(vGridCapital, tCapitalPolicyAccelerator[:,end,end], label = "z_5, A_3")
plot!(vGridCapital, vGridCapital, color=:black, linestyle=:dash, label= "45 degrees")

savefig(pCapitalPolicyAccelerator, "/Users/Castesil/Documents/EUI/Year II - PENN/Computational Economics/Homework 1/LaTeX/pCapitalPolicyAccelerator.pdf")

# 6.2
pConvergenceAccelerator = plot(vMaxDifferenceAccelerator[2:end], ylabel = "error", xlabel = "number of iterations", legend=:none)
savefig(pConvergenceAccelerator, "/Users/Castesil/Documents/EUI/Year II - PENN/Computational Economics/Homework 1/LaTeX/pConvergenceAccelerator.pdf")

# --------------
# 7. Multigrid
# --------------
nMidPoints = [7, 9, 10]  # construct more grid so built up dict
# approx nGridCapital equal 7=100, 9=500, 12=5000

@time tValueFunction, tCapitalPolicy, vGridCapital, vMaxDifference = exercise7.multi_grid(villaverdeEconomy, villaverdeSteadyState, nMidPoints)

pCapitalPolicyMultigrid = plot(vGridCapital, tCapitalPolicy[:,1,1], label = "z_1, A_1", xlabel = "capital", legend=:bottomright)
plot!(vGridCapital, tCapitalPolicy[:,end,1],label = "z_5, A_1")
plot!(vGridCapital, tCapitalPolicy[:,1,end], label = "z_1, A_5")
plot!(vGridCapital, tCapitalPolicy[:,end,end], label = "z_5, A_3")
plot!(vGridCapital, vGridCapital, color=:black, linestyle=:dash, label= "45 degrees")

savefig(pCapitalPolicyMultigrid, "/Users/Castesil/Documents/EUI/Year II - PENN/Computational Economics/Homework 1/LaTeX/pCapitalPolicyMultigrid.pdf")

# -------------------
# 8. Stochastic grid
# -------------------
