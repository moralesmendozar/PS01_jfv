module ex3c
export c_multigrid_enhanced
using Parameters, LinearAlgebra, Interpolations, JLD2, NLsolve

include("gridNum.jl")
using .ex7aux

function c_multigrid_enhanced(economy, steadyStateValues, nMidPoints::Array)
    # 0. Get (unpack) the parameters back
    @unpack vGridZ, vGridA,mTranstnZ,mTranstnA, mTranstnZA,α,β,δ,θ,maxiter = economy
    @unpack kss, l1ss, l2ss, utilitySS = steadyStateValues
    nZ = length(vGridZ)
    nA = length(vGridA)

    # Grid Capital
    kmin = 0.7 * kss
    kmax = 1.3 * kss

    dicLabor = Dict()
    dictVF = Dict()

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

    # --------------------------------------------------------------------------
    # Multigrid scheme calculation... :)
    println(" Multigrid starts ;).... ")
    tol = 1.0e-6
    global iteration = 0

    for iGrid in 1:length(nMidPoints)

        println("Solving for grid number $iGrid ")
        println(" ")
        # -----------------
        # Pre-allocation for VFI stuff
        if iGrid == 1
            # Generate the Capital grid:
            global vGridK = ex7aux.gridNum(kmin, kmax, nMidPoints[iGrid])
            nK = length(vGridK)
            global tVF = fill(utilitySS, nK, nZ, nA)
        else
            itp = interpolate((vGridK, vGridZ, vGridA), tVF, Gridded(Linear()))
            global vGridK = ex7aux.gridNum(kmin, kmax, nMidPoints[iGrid])
            tVF = itp(vGridK, vGridZ,vGridA)
        end
        nK = length(vGridK)
        tVFNew = zeros(nK, nZ, nA)
        global tPolicyFnIndx = zeros(Int64,nK,nZ, nA)
        global d4ProvisionalLaborOne = zeros(nK,nK, nZ, nA)
        global d4ProvisionalLaborTwo = zeros(nK,nK, nZ, nA)
        global vProvisional = zeros(nK)
        global vMaxDiff = Float64[]
        # ----------------------------------------------------------------------
        # Find Labour choice matrices
        # Initial guess for labour choices
        guessL1 = l1ss
        guessL2 = l2ss

        # IMPROVE THE MATRIX FOR LABOR CALCULATION, ENHANCE...
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
                        keyLabor = (vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA])
                        if keyLabor ∈ keys(dicLabor)
                            l1, l2 = dicLabor[keyLabor]

                        elseif l1prev == l1ss
                            l1 = l1ss
                            l2 = l2ss
                            keyLabor = (vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA])
                            dicLabor[keyLabor] = [l1, l2]

                        else

                            l1, l2 = try
                                labour_choice(vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA], guessL1, guessL2,α,β,δ,θ,eZ)
                            catch
                                guessL1, guessL2
                            end

                            # store in dictionary
                            keyLabor = (vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA])
                            dicLabor[keyLabor] = [l1, l2]
                        end #catch key_labor

                        l1prev = l1
                        d4ProvisionalLaborOne[iK,iKNext,iZ,iA] = l1
                        d4ProvisionalLaborTwo[iK,iKNext,iZ,iA] = l2

                    end #for iKNext in iK:nK
                    #do left triangle
                    if iK>=2
                        for iKNext in iK-1:-1:1

                            # Optimal labor choices
                            keyLabor = (vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA])
                            if keyLabor ∈ keys(dicLabor)
                                l1, l2 = dicLabor[keyLabor]

                            elseif l2prev == l2ss
                                l1 = l1ss
                                l2 = l2ss
                                keyLabor = (vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA])
                                dicLabor[keyLabor] = [l1, l2]

                            else

                                l1, l2 = try
                                    labour_choice(vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA], guessL1, guessL2,α,β,δ,θ,eZ)
                                catch
                                    guessL1, guessL2
                                end

                                # store in dictionary
                                keyLabor = (vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA])
                                dicLabor[keyLabor] = [l1, l2]
                            end #catch key_labor

                            l2prev = l2
                            d4ProvisionalLaborOne[iK,iKNext,iZ,iA] = l1
                            d4ProvisionalLaborTwo[iK,iKNext,iZ,iA] = l2

                        end #for iKNext in iK:nK
                    end
                end  #for iK in 1:nK
            end  #for iZ in 1:nZ
        end  #for iA in 1:nA

        # ----------------------------------------------------------------------
        # VFI (computed each multigridIteration)
        global maxDiff = 1.0
        iteration = 0

        println("VFI ...")
        println(" ")

        while(maxDiff > tol)
            for iA in 1:nA
                for iZ in 1:nZ
                    for iK in 1:nK
                        if iteration ==0
                            #check keyVF
                            keyValueFunction = (vGridK[iK], vGridZ[iZ], vGridA[iA])
                            if keyValueFunction ∈ keys(dictVF)
                                tVFNew[iK, iZ, iA], tPolicyFnIndx[iK, iZ, iA] = dictVF[keyValueFunction]
                            else
                                #compute Vprovisional
                                for iKNext in 1:nK

                                    expected = 0.0

                                    for iAnext in 1:nA
                                        for iZnext in 1:nZ
                                            expected = expected + mTranstnZ[iZ, iZnext] * mTranstnA[iA, iAnext] * tVF[iKNext,iZnext,iAnext]
                                        end
                                    end

                                    # labor
                                    keyLabor = (vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA])
                                    l1, l2 = dicLabor[keyLabor]
                                    # Values for all k'
                                    vProvisional[iKNext] = val_provisional(vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA], l1, l2, expected, economy)

                                end #for iKNext in 1:nK

                                tVFNew[iK, iZ, iA], tPolicyFnIndx[iK, iZ, iA] = findmax(vProvisional)
                                dictVF[keyValueFunction] = [tVFNew[iK, iZ, iA], tPolicyFnIndx[iK, iZ, iA]]
                            end #if check keyVF
                        else
                            #compute valProvisional
                            for iKNext in 1:nK

                                expected = 0.0

                                for iAnext in 1:nA
                                    for iZnext in 1:nZ
                                        expected = expected + mTranstnZ[iZ, iZnext] * mTranstnA[iA, iAnext] * tVF[iKNext,iZnext,iAnext]
                                    end
                                end

                                # labor
                                keyLabor = (vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA])
                                l1, l2 = dicLabor[keyLabor]
                                # Values for all k'
                                vProvisional[iKNext] = val_provisional(vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA], l1, l2, expected, economy)

                            end #for iKNext in 1:nK

                            tVFNew[iK, iZ, iA], tPolicyFnIndx[iK, iZ, iA] = findmax(vProvisional)

                            #... and save the VF in the dictionary for VFI
                            keyValueFunction = (vGridK[iK], vGridZ[iZ], vGridA[iA])
                            dictVF[keyValueFunction] = [tVFNew[iK, iZ, iA], tPolicyFnIndx[iK, iZ, iA]]

                        end #check iteration ==1

                    end #for iK in 1:nK
                end #for iZ in 1:nZ
            end #for iA in 1:nA

            maxDiff = norm(tVFNew-tVF)
            vMaxDiff = push!(vMaxDiff,maxDiff)

            tVF, tVFNew = tVFNew, tVF

            iteration = iteration+1
            # if(mod(iteration,10)==0 || iteration == 1)
            #     println(" Iteration = ", iteration, " Sup Diff = ", maxDiff)
            # end #if
        end #while VFI

        println("Multigrid for grid number $iGrid finished solving")
        println(" Number of Iterations = ", iteration, " Sup Diff = ", maxDiff)
    end  #for multigrid...

    println("Multigrid finished")
    println(" Number of Iterations = ", iteration, " Sup Diff = ", maxDiff)

    tPolicyFn = vGridK[tPolicyFnIndx]
    #tPolicyFn = tPolicyFnIndx
    return tVF, tPolicyFn, vGridK, vMaxDiff, d4ProvisionalLaborOne, d4ProvisionalLaborTwo
end

end
