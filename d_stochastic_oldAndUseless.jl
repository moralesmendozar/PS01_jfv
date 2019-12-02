module ex3d
export d_stochastic
using Parameters, LinearAlgebra, Interpolations, JLD2, NLsolve, Random

include("gridNum.jl")
using .ex7aux

function d_stochastic(economy, steadyStateValues,vGridKComplete, nKComplete,Vinit)
    # 0. Get (unpack) the parameters back
    @unpack vGridZ, vGridA,mTranstnZ,mTranstnA, mTranstnZA,α,β,δ,θ,maxiter = economy
    @unpack kss, l1ss, l2ss, utilitySS = steadyStateValues
    nZ = length(vGridZ)
    nA = length(vGridA)

    # Grid Capital
    kmin = vGridKComplete[1]
    kmax = vGridKComplete[end]

    # Auxiliary function
    function labour_choice(K,Knext,ProdZ, ProdA, guessL1, guessL2,α,β,δ,θ,eZ)
        # FUNCTION FINDS THE L1,L2 GIVEN K,K', using the FOCs of SPP.... :)
        #       using the nlsolve package (function)  :)

        function focs12!(F,x)
            c1 = eZ * K^α * x[1]^(1-α) + (1-δ)*K - Knext
            # FOC [l1]
            F[1] = θ * (1-α)*eZ*K^α*x[1]^(-α) * (c1)^(θ-1)*(ProdA*x[2])^(1-θ) - x[1] - x[2]
            # FOC [l2]
            F[2] = ( 1 - θ)* (c1)^θ *(ProdA*x[2])^(-θ)*ProdA - x[1] - x[2]
        end

        solvels = nlsolve(focs12!, [guessL1, guessL2])
        LOne = solvels.zero[1]
        LTwo = solvels.zero[2]

        return LOne, LTwo
    end

    function getLabourOrReturnGuess(gridK, gridKNext, gridProdZ, gridProdA, guess1, guess2,α,β,δ,θ,eZ)
        laborOne, laborTwo = try
            labour_choice(gridK, gridKNext, gridProdZ, gridProdA, guess1, guess2,α,β,δ,θ,eZ)
        catch  #make labor = 0 if there is no solution
            guess1, guess2
        end
        return laborOne, laborTwo
    end

    function val_provisional(gridK::Real, gridKNext::Real, gridProdZ::Real, gridProdA::Real, laborOne::Real, laborTwo::Real, eV::Real, econ)
        # Function gets l1,l2
        #       (and the continuation value given k,k' l1(k,k'),l2(k,k'))

        # Get parameters:
        @unpack α, β, δ, θ = econ
        eZ = exp(gridProdZ)

        # Get consumptions given the labours
        c1 = eZ*gridK^α * laborOne^(1-α) + (1-δ) * gridK - gridKNext
        c2 = gridProdA * laborTwo

        # Continuation value
        if c1 < 0 || c2 < 0
            valueProvisional = -100000
        else
            ut = c1^θ*c2^(1-θ) - 0.5*(laborOne+laborTwo)^2
            valueProvisional = (1-β)*ut + β*eV
        end

        return valueProvisional
    end


    # Find Labour choice matrices
    global mL1 = zeros(nK,nK, nZ, nA)
    global mL2 = zeros(nK,nK, nZ, nA)
    # Initial guess for labour choices
    guessL1 = l1ss
    guessL2 = l2ss

    # MATRIX FOR LABOR CALCULATION, ENHANCED...
    println("Finding labour choice matrices for all k,kprime' ... ")
    println(" ")
    for iA in 1:nA
        for iZ in 1:nZ
            #reset l1prev and l2prev
            l1prev = -10
            l2prev = -10
            for iK in 1:nK
                #do right triangle
                for iKNext in iK:nK

                    # Optimal labor choices
                    if l1prev == l1ss
                        l1 = l1ss
                        l2 = l2ss
                    else
                        l1, l2 = try
                            labour_choice(vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA], guessL1, guessL2,α,β,δ,θ,eZ)
                        catch
                            guessL1, guessL2
                        end
                    end
                    l1prev = l1
                    mL1[iK,iKNext,iZ,iA] = l1
                    mL2[iK,iKNext,iZ,iA] = l2

                end #for iKNext in iK:nK
                #do left triangle
                if iK>=2
                    for iKNext in iK-1:-1:1

                        # Optimal labor choices
                        if l2prev == l2ss
                            l1 = l1ss
                            l2 = l2ss
                        else
                            l1, l2 = try
                                labour_choice(vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA], guessL1, guessL2,α,β,δ,θ,eZ)
                            catch
                                guessL1, guessL2
                            end
                        end #catch key_labor
                        l2prev = l2
                        mL1[iK,iKNext,iZ,iA] = l1
                        mL2[iK,iKNext,iZ,iA] = l2

                    end #for iKNext in iK:nK
                end
            end  #for iK in 1:nK
        end  #for iZ in 1:nZ
    end  #for iA in 1:nA

    #make dictionary mapping iik to [iK,iZ,iA]
    iik = 0
    dicIkzaReversed = Dict()
    vecIkza = zeros(nK*nZ*nA,3)
    for iA in 1:nA
        for iZ in 1:nZ
            for iK in 1:nK
                iik = iik +1
                vecIkza[iik,:] = [iK,iZ,iA]
                #make dictionary with the reverse:
                keyIk = [iK,iZ,iA]
                dicIkzaReversed[keyIk] = iik
            end  #for iK in 1:nK
        end  #for iZ in 1:nZ
    end  #for iA in 1:nA


    #STOCHASTIC GRID SCHEME....
    println(" Stochastic grid starts ;).... ")
    tol = 1.0e-6
    kiter = 0
    maxdif = 10.0
    N0 = 10
    #dictVF = Dict() #unsure if will be needed....
    tVF = Vinit
    global ip = interpolate((vGridKComplete, vGridZ, vGridA), tVF, Gridded(Linear()))
    global policyfnip = ip
    global vMaxDiff = Float64[]

    while maxdif > tol
        Nkiter = 2^(2*kiter)*N0
        kiter = kiter+1

        #gets a sample from all possible points
        Spk = Random.Sampler(MersenneTwister(kiter), 1:(nKComplete*nA*nZ))
        SkI = rand(rng, sp,Nkiter,1)  #gets Nkiter-size sample from the possible points
        SkI = sort(SkI)
        sis = vecIkza[SkI] #contains all the points to explore, in [ik,iz,iA] form
        #       some could be repeated
        #Find the probability matrix P([iZ,iA],ik,->[iz',ia'])
        PP = zeros(nZ,nA,nK,nKiter)
        for iA in 1:nA
            for iZ in 1:nZ
                for ikp in 1:nK
                    suma = 0
                    for ii in SkI
                        kii = sis[ii,1]
                        zii = sis[ii,2]
                        Aii = sis[ii,3]
                        if ik == kii
                            PP[iZ,iA,ikp,ii] = mTranstnZ[iZ, zii] * mTranstnA[iA, Aii]
                            suma = suma + PP[iZ,iA,ik,ikp,ii]
                        end
                    end
                    PP[iZ,iA,ikp,:] = PP[iZ,iA,ikp,:]./suma
                end
            end
        end #loop for PP probabilities

        ############################ check interpolations :O
        VPrevious = ip(sis[:,1],sis[:,2],sis[:,3])
        #update, find max
        VFhat, policyK = gamma(VPrevious)










        maxDiff = norm(VFhat-VPrevious)
        vMaxDiff = push!(vMaxDiff,maxDiff)
        ############################ check interpolations :O
        #interpolation is problem #1
        ip = interpolate((sis[:,1],sis[:,2],sis[:,3]),Vfhat)
        policyfnip = interpolate(sis[:,1],sis[:,2],sis[:,3],policyK)

    end #while dif>tol

    #tPolicyFn = vGridK[tPolicyFnIndx]
    #tPolicyFn = tPolicyFnIndx
    return ip, policyfnip, vGridKComplete, vMaxDiff, mL1, mL2
end

end
