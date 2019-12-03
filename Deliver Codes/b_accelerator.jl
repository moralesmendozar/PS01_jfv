module ex03b
export b_accelerator
using Parameters, LinearAlgebra, Interpolations, JLD2

function b_accelerator(economy, steadyStateValues,Vinit,nK)
    # 0. Get (unpack) the parameters back
    @unpack vGridZ, vGridA,mTranstnZ,mTranstnA, mTranstnZA,α,β,δ,θ,maxiter = economy
    @unpack kss, l1ss, l2ss = steadyStateValues
    nZ = length(vGridZ)
    nA = length(vGridA)

    # Grid Capital
    vGridK = collect(range(0.7 * kss, 1.3 * kss, length = nK))

    # Pre-allocation
    #   Create empty matrices where to record results (to return also) (:)
    tVF = Vinit
    tVFNew = zeros(nK,nZ,nA)
    tPolicyFnIndx = zeros(Int64,nK,nZ,nA)
    vProvisional = zeros(nK)

    # Initial guess for labour choices
    guessL1 = l1ss
    guessL2 = l2ss

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

    # Get utility(z,a,k) and labor(z,A,k,k') functions (matrices, actually)
    #    for all (z,A,k,k')
    #       to use in VFI to find V(z,A,k) and k'(z,A,k)
    println(" Getting tensor of labors... ")
    ml1 = -22*ones(nK,nK,nZ,nA)
    ml2 = -22*ones(nK,nK,nZ,nA)
    # for iA in 1:nA
    #     for iZ in 1:nZ
    #         for iK in 1:nK
    #             for iKNext in 1:nK
    #                 eZ = exp(vGridZ[iZ])
    #                 l1sol, l2sol = getLabourOrReturnGuess(vGridK[iK], vGridK[iKNext], vGridZ[iZ], vGridA[iA], l1ss, l2ss,α,β,δ,θ,eZ)
    #                 ml1[iK,iKNext,iZ,iA] = l1sol
    #                 ml2[iK,iKNext,iZ,iA] = l2sol
    #             end
    #         end
    #     end
    # end

    for iA in 1:nA
        for iZ in 1:nZ
            eZ = exp(vGridZ[iZ])
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
                    ml1[iK,iKNext,iZ,iA] = l1
                    ml2[iK,iKNext,iZ,iA] = l2

                end #for iKNext in iK:nKComplete
                #do left triangle
                if iK>=2
                    for iKNext in iK-1:-1:1

                        # Optimal labor choices
                        if l2prev == l2ss
                            l1 = l1ss
                            l2 = l2ss
                        else
                            l1, l2 = try
                                labour_choice(vGridKComplete[iK], vGridKComplete[iKNext], vGridZ[iZ],vGridA[iA], guessL1, guessL2,α,β,δ,θ,eZ)
                            catch
                                guessL1, guessL2
                            end
                        end #catch key_labor
                        l2prev = l2
                        ml1[iK,iKNext,iZ,iA] = l1
                        ml2[iK,iKNext,iZ,iA] = l2

                    end #for iKNext in iK:nKComplete
                end
            end  #for iK in 1:nKComplete
        end  #for iZ in 1:nZ
    end  #for iA in 1:nA


    println(" Tensor of labors computed... ")

    println(" VFI starts.... ")
    # VFI
    maxDiff = 10.0
    tol = 1.0e-6
    global iteration = 0
    global vMaxDiff = Float64[]

    while(maxDiff > tol && iteration<maxiter)

        #itp = interpolate((vGridK, vGridZ, vGridA), tVF, Gridded(Linear()))
        #tVFInterpol = itp(vGridK, vGridZ,vGridA)

        #ACCELERATOR, do max only once every ten
        if(mod(iteration,10)==0 || iteration == 1)
            #do max step only once every ten times:
            for iA in 1:nA
                for iZ in 1:nZ
                    for iK in 1:nK
                        for iKNext in 1:nK
                            #get expected value:
                            expected = 0.0
                            for iAnext in 1:nA
                                for iZnext in 1:nZ
                                    expected = expected + mTranstnZ[iZ, iZnext] * mTranstnA[iA, iAnext] * tVF[iKNext,iZnext,iAnext]
                                end
                            end #for iAnext

                            eZ = exp(vGridZ[iZ])
                            L1solved = ml1[iK,iKNext,iZ,iA]#vGridK[iK]^(1-α)
                            L2solved = ml2[iK,iKNext,iZ,iA]#vGridK[iK]^(1-α)
                            vProvisional[iKNext] = val_provisional(vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA],L1solved, L2solved, expected, economy)
                        end #for iKNext

                        tVFNew[iK, iZ, iA], tPolicyFnIndx[iK, iZ, iA] = findmax(vProvisional)
                    end #for iK
                end # for iZ
            end  #end of for iA and fors...

        else
            #Do accelerator step, just update VF
            for iA in 1:nA
                for iZ in 1:nZ
                    for iK in 1:nK


                        iKnextOptimal = tPolicyFnIndx[iK, iZ, iA] #found previously
                        kNxt = vGridK[iKnextOptimal]
                        #get expected value:
                        expectedVal = 0.0
                        for iAnext in 1:nA
                            for iZnext in 1:nZ
                                expectedVal = expectedVal + mTranstnZ[iZ, iZnext] * mTranstnA[iA, iAnext] * tVF[iKnextOptimal,iZnext,iAnext]
                            end
                        end #for iAnext

                        l1 = ml1[iK,iKnextOptimal,iZ,iA]#vGridK[iK]^(1-α)
                        l2 = ml2[iK,iKnextOptimal,iZ,iA]#vGridK[iK]^(1-α)

                        eZ = exp(vGridZ[iZ])
                        kNow = vGridK[iK]
                        c1 = eZ*kNow^α * l1^(1-α) + (1-δ) * kNow - kNxt
                        c2 = vGridA[iA] * l2
                        if c1 < 0 || c2 < 0
                            valueProvisional = -100000
                        else
                            util = c1^θ*c2^(1-θ) - 0.5*(l1+l2)^2
                            valueProvisionalV = (1-β)*util + β*expectedVal
                        end

                        tVFNew[iK, iZ, iA],  = valueProvisionalV

                    end #for iK
                end # for iZ
            end  #end of for iA and fors...
        end #end of if (for accelerator comparison)

        maxDiff = norm(tVFNew-tVF)
        vMaxDiff = push!(vMaxDiff,maxDiff)
        tVF, tVFNew = tVFNew, tVF

        iteration = iteration+1
        # if(mod(iteration,10)==0 || iteration == 1)
        #     println(" Iteration = ", iteration, " Sup Diff = ", maxDiff)
        # end #end of if for the println(iteraton and sup diff)
    end  #end of while()
    println("Accelerator finished")
    println(" Number of Iterations = ", iteration, " Sup Diff = ", maxDiff)
    tPolicyFn = vGridK[tPolicyFnIndx]
    return tVF, tPolicyFn, vGridK, vMaxDiff
end

end
