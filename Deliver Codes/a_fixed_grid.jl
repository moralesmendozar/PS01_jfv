module ex3
export a_fixed_grid
using Parameters, LinearAlgebra, Interpolations

function a_fixed_grid(economy, steadyStateValues,Vinit,nK)

    # 0. Get (unpack) the parameters back
    @unpack vGridZ, vGridA,mTranstnZ,mTranstnA, mTranstnZA,α,β,δ,θ,maxiter = economy
    @unpack kss, l1ss, l2ss = steadyStateValues
    nZ = length(vGridZ)
    nA = length(vGridA)

    # Grid Capital
    vGridK = collect(range(0.7 * kss, 1.3 * kss, length = nK))

    # Pre-allocation
    #   Create empty matrices where to record results (to return also) (:)
    #tVF = zeros(nK,nZ,nA)
    #tVF = repeat(vGridK, nZ,nA)
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
        #function parameters!(F,x, K,Knext,α,δ,eZ)
        #    c1 = eZ * K^α * x[1]^(1-α) + (1-δ)*K - Knext
        #end

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

    function val_provisional(gridK::Real, gridKNext::Real, gridProdZ::Real, gridProdA::Real, guess1::Real, guess2::Real, eV::Real, econ)
        # Function gets l1,l2
        #       (and the continuation value given k,k' l1(k,k'),l2(k,k'))

        # Get parameters:
        @unpack α, β, δ, θ = econ
        eZ = exp(gridProdZ)

        # Find Labor l1,l2 choices, found given K,K':
        laborOne, laborTwo = try
            labour_choice(gridK, gridKNext, gridProdZ, gridProdA, guess1, guess2,α,β,δ,θ,eZ)
        catch  #make labor = 0 if there is no solution
            guess1, guess2
        end

        # Get consumptions given the found labours
        c1 = eZ*gridK^α * laborOne^(1-α) + (1-δ) * gridK - gridKNext
        c2 = gridProdA * laborTwo

        # 3. Value
        if c1 < 0 || c2 < 0
            valueProvisional = -100000
        else
            ut = c1^θ*c2^(1-θ) - 0.5*(laborOne+laborTwo)^2
            valueProvisional = (1-β)*ut + β*eV
        end

        return valueProvisional, laborOne, laborTwo
    end

    # Get utility(z,a,k) and labor(z,A,k,k') functions (matrices, actually)
    #    for all (z,A,k,k')
    #       to use in VFI to find V(z,A,k) and k'(z,A,k)
    ml1 = zeros(nK,nZ,nA)



    # VFI
    maxDiff = 10.0
    tol = 1.0e-6
    iteration = 0

    while(maxDiff > tol && iteration<maxiter)

        itp = interpolate((vGridK, vGridZ, vGridA), tVF, Gridded(Linear()))
        tVFInterpol = itp(vGridK, vGridZ,vGridA)

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
                        end

                        eZ = exp(vGridZ[iZ])
                        #vProvisionalValues[iCapitalNext] = val_provisional(vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA], expected, economy)
                        guessL1 = l1ss#vGridK[iK]^(1-α)
                        guessL2 = l2ss#vGridK[iK]^(1-α)
                        vProvisional[iKNext], l1, l2 = val_provisional(vGridK[iK], vGridK[iKNext], vGridZ[iZ],vGridA[iA],guessL1, guessL2, expected, economy)
                        #guessL1 = l1
                        #guessL2 = l2

                    end
                    tVFNew[iK, iZ, iA], tPolicyFnIndx[iK, iZ, iA] = findmax(vProvisional)
                end
            end
        end  #end of fors...

        maxDiff = norm(tVFNew-tVF)
        tVF, tVFNew = tVFNew, tVF

        iteration = iteration+1
        if(mod(iteration,10)==0 || iteration == 1)
            println(" Iteration = ", iteration, " Sup Diff = ", maxDiff)
        end #end of if
    end  #end of while()

    tPolicyFn = vGridK[tPolicyFnIndx]
    return tVF, tPolicyFn, vGridK

end

end
