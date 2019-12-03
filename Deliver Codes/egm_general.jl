module ex4

export egm_general

# Packages
using Parameters
using LinearAlgebra, Interpolations, NLsolve

# Modules
include("./egm_lss.jl")
using .ex4ss
include("./egm_vfi_std.jl")
using .ex4vfi

function egm_general(economy, steadyState, nK)
    # Parameters
    @unpack α, θ, δ, vGridZ, mTranstnZ, vGridA, mTranstnA = economy
    @unpack kss, l1ss, l2ss, utilitySS = steadyState

    nGridZ = length(vGridZ)
    nGridA = length(vGridA)

    # Kapital Grid
    vGridK = collect(range(0.7 * kss, 1.3 * kss, length = nK))

    # 1. Guess VFI (initial)
    # make it increasing
    vInitialGuess = Float64[]
    for iK in 1:nK
        value = (iK-round(nK/2))/nK + utilitySS
        push!(vInitialGuess,value)
    end

    # 2. EGM: holding labor in its steady state value
    tValueFunctionTilde, tMarketResourcesEndogenous = ex4ss.egm_lss(economy, steadyState, vGridK, vInitialGuess)

    # 3. Recover opitmal policies
    #tPolicyFunction, tConsPolicy = ex4vfi.egm_vfi_std(tValueFunctionTilde, tMarketResourcesEndogenous, economy, steadyState, vGridK)
    tVFtilde, tPolicyFunction, tCons1Policy, tCons2Policy, tL1Policy, tL2Policy = ex4vfi.egm_vfi_std(tValueFunctionTilde,  tMarketResourcesEndogenous, economy, steadyState, vGridK)
    #return tValueFunctionTilde, tPolicyFunction, tConsPolicy
    return tValueFunctionTilde, tVFtilde, tPolicyFunction, tCons1Policy, tCons2Policy, tL1Policy, tL2Policy
end

end
