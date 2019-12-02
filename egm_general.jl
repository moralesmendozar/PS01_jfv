module exercise4

export egm_labor

# ---------
# Packages
# ---------
using Parameters
using LinearAlgebra, Interpolations, NLsolve

# --------
# Modules
# --------
include("./egm_lss.jl")
using .exercise4steadystate

include("./egm_vfi_standard.jl")
using .exercise4vfi

function egm_labor(economy, steadyState, nGridCapital)

    # ---------------
    # 0. Housekeeping
    # ---------------
    @unpack α, θ, δ, vGridZ, mTransitionZ, vGridA, mTransitionA = economy
    @unpack capitalSteadyState, laborOneSteadyState, laborTwoSteadyState, utilitySteadyState = steadyState

    nGridZ = length(vGridZ)
    nGridA = length(vGridA)

    # ----------------
    # 1. Grid Capital
    # ----------------
    vGridCapitalNext = collect(range(0.7 * capitalSteadyState, 1.3 * capitalSteadyState, length = nGridCapital))

    # --------------------------------
    # 2. Guess for the value function
    # --------------------------------
    # needs to be increasing
    vInitialGuess = Float64[]
    for iCapital in 1:nGridCapital
        value = iCapital/nGridCapital + utilitySteadyState
        push!(vInitialGuess,value)
    end

    # -------------------------------------------------
    # 3. EGM: holding labor to its steady state value
    # -------------------------------------------------
    tValueFunctionTilde = exercise4steadystate.egm_labor_ss(economy, steadyState, vGridCapitalNext, vInitialGuess)

    # ----------------------------
    # 4. Recover opitmal policies
    # -----------------------------
    tValueFunctionTilde = exercise4vfi.egm_vfi_standard(economy, steadyState, vGridCapitalNext, vInitialGuess)
