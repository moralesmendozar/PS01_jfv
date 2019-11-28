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

# ---------------
# 1. Calibration
# ---------------
villaverde_economy = @with_kw (
                        # a) Preferences
                        β = 0.96, # time discount factor
                        θ = 0.5,  # elasticity of substitution between goods

                        # b) Technology
                        α = 0.33, # capital share
                        δ = 0.9, # depreciation

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
@unpack α, δ = villaverdeEconomy

# Solving system of equations
results = exercise2.ss_solve()
capitalSteadyState = results.zero[1]
laborOneSteadyState = results.zero[2]
laborTwoSteadyState = results.zero[3]

# From the resource constraints
consumptionOneSteadyState = capitalSteadyState^(α) * laborOneSteadyState^(1-α) - (1 - δ) * capitalSteadyState
consumptionTwoSteadyState = laborTwoSteadyState


villaverde_SteadyState = @with_kw (
                        capitalSteadyState = capitalSteadyState,

                        laborOneSteadyState = laborOneSteadyState,

                        consumptionOneSteadyState = consumptionOneSteadyState,

                        laborTwoSteadyState = laborTwoSteadyState,

                        consumptionTwoSteadyState = consumptionTwoSteadyState,

                        utilitySteadyState = consumptionOneSteadyState^(0.5) * consumptionTwoSteadyState^(0.5) - 0.5*(laborOneSteadyState + laborTwoSteadyState)^2
)

villaverdeSteadyState = villaverde_SteadyState()

println(" Deterministic steady state ")
println(" -------------------------- ")
println(" K_ss = ", round(capitalSteadyState, digits= 5), "  c_ss^1 = ",round(consumptionOneSteadyState, digits=5), "  l_ss^1 = ", round(laborOneSteadyState, digits=5), "  c_ss^2 = ", round(consumptionTwoSteadyState, digits=5), "  l_ss^2 = ", round(laborTwoSteadyState, digits=5))
println(" ")

# --------------
# 3. Fixed grid
# --------------

@time tValueFunction, tPolicyFunction, vGridCapital = exercise3.fixed_grid(villaverdeEconomy, villaverdeSteadyState)

# Graphs
pValueFunction = plot(vGridCapital, tValueFunction[:,1,1], label = "z_1, A_1", xlabel = "Capital", legend=:topright)
plot!(vGridCapital, tValueFunction[:,end,1], label = "z_5, A_1")
plot!(vGridCapital, tValueFunction[:,1,end], label = "z_1, A_5")
plot!(vGridCapital, tValueFunction[:,end,end], label = "z_5, A_3")


pPolicyFunction = plot(vGridCapital, tPolicyFunction[:,1,1], label = "z_1, A_1", xlabel = "Capital", color=:blue)
plot!(vGridCapital, tPolicyFunction[:,end,1],label = "z_5, A_1", color=:blue, linestyle=:dash)
plot!(vGridCapital, tPolicyFunction[:,1,end], label = "z_1, A_5", color=:red)
plot!(vGridCapital, tPolicyFunction[:,end,end], color=:red, linestyle=:dash, label = "z_5, A_3")
plot!(vGridCapital,vGridCapital, color=:black,linestyle=:dash)
