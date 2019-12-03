module exercise4

export egm_labor

# ---------
# Packages
# ---------
using Parameters
using LinearAlgebra, Interpolations, NLsolve, Optim

# --------
# Modules
# --------
include("./src/egm_lss.jl")
using .exercise4steadystate

include("./src/egm_vfi_standard.jl")
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
        value = (iCapital-(nGridCapital÷2))/nGridCapital + utilitySteadyState
        push!(vInitialGuess,value)
    end

    # -------------------------------------------------
    # 3. EGM: holding labor to its steady state value
    # -------------------------------------------------
    tValueFunctionTilde, tMarketResourcesEndogenous = exercise4steadystate.egm_labor_ss(economy, steadyState, vGridCapitalNext, vInitialGuess)

    # ----------------------------
    # 4. Recover opitmal policies
    # -----------------------------
    tPolicyFunction, tConsumptionOnePolicy, tConsumptionTwoPolicy, tLaborOnePolicy, tLaborTwoPolicy = exercise4vfi.egm_vfi_standard(tValueFunctionTilde,  tMarketResourcesEndogenous, economy, steadyState, vGridCapitalNext)

    # ------------------------------------------------
    # 5. Apply the EGM using new labor policy function
    # ------------------------------------------------

    # 5.1 Find values of k s.t. k' = g(k,Z,A) = nGridCapitalNext
    #while maxDifference > tolerance


    return tValueFunctionTilde, vGridCapitalNext, tPolicyFunction, tConsumptionOnePolicy, tConsumptionTwoPolicy, tLaborOnePolicy, tLaborTwoPolicy

end

end
