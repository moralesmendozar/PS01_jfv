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
@unpack α, β, δ = villaverdeEconomy

capitalSteadyState =  ((α * β)/(1 - β * δ))^(1/(1-α))

laborSteadyState = (1/1-α)*(1-(1-δ)*capitalSteadyState^(1-α)) + 1

laborTwoSteadyState = laborSteadyState - 1

consumptionOneSteadyState = capitalSteadyState^α - (1-δ)*capitalSteadyState

consumptionTwoSteadyState = laborTwoSteadyState

consumptionSteadyState = consumptionOneSteadyState + consumptionTwoSteadyState

println(" Deterministic steady state ")
println(" -------------------------- ")
println(" K_ss = ", capitalSteadyState, "  L_ss = ", laborSteadyState, "  C_ss = ", consumptionSteadyState)
println(" ")

# --------------
# 3. Fixed grid
# --------------

@time tValueFunction, tPolicyFunction, vGridCapital = exercise3.fixed_grid(villaverdeEconomy, capitalSteadyState)

# Graphs
pValueFunction = plot(vGridCapital, tValueFunction[:,1,1], label = "z_1, A_1", xlabel = "Capital", legend=:topleft)
plot!(vGridCapital, tValueFunction[:,end,1], label = "z_5, A_1")
plot!(vGridCapital, tValueFunction[:,1,end], label = "z_1, A_5")
plot!(vGridCapital, tValueFunction[:,end,end], label = "z_5, A_3")


pPolicyFunction =  plot(vGridCapital,vGridCapital, color=:black,linestyle=:dash)
plot!(vGridCapital, tPolicyFunction[:,1,1], label = "z_1, A_1", xlabel = "Capital", color=:blue)
plot!(vGridCapital, tPolicyFunction[:,end,1],label = "z_5, A_1", color=:blue, linestyle=:dash)
plot!(vGridCapital, tPolicyFunction[:,1,end], label = "z_1, A_5", color=:red)
plot!(vGridCapital, tPolicyFunction[:,end,end], color=:red, linestyle=:dash, label = "z_5, A_3")
