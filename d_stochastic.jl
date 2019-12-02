module ex3d
export d_stochastic
using Parameters, LinearAlgebra, Interpolations, JLD2, NLsolve, Random

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
    global mL1 = zeros(nKComplete,nKComplete, nZ, nA)
    global mL2 = zeros(nKComplete,nKComplete, nZ, nA)
    # Initial guess for labour choices
    guessL1 = l1ss
    guessL2 = l2ss

    # MATRIX FOR LABOR CALCULATION, ENHANCED...
    println("Finding labour choice matrices for all k,kprime' ... ")
    println(" ")
    for iA in 1:nA
        for iZ in 1:nZ
            eZ = exp(vGridZ[iZ])
            #reset l1prev and l2prev
            l1prev = -10
            l2prev = -10
            for iK in 1:nKComplete
                #do right triangle
                for iKNext in iK:nKComplete

                    # Optimal labor choices
                    if l1prev == l1ss
                        l1 = l1ss
                        l2 = l2ss
                    else
                        l1, l2 = try
                            labour_choice(vGridKComplete[iK], vGridKComplete[iKNext], vGridZ[iZ],vGridA[iA], guessL1, guessL2,α,β,δ,θ,eZ)
                        catch
                            guessL1, guessL2
                        end
                    end
                    l1prev = l1
                    mL1[iK,iKNext,iZ,iA] = l1
                    mL2[iK,iKNext,iZ,iA] = l2

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
                        mL1[iK,iKNext,iZ,iA] = l1
                        mL2[iK,iKNext,iZ,iA] = l2

                    end #for iKNext in iK:nKComplete
                end
            end  #for iK in 1:nKComplete
        end  #for iZ in 1:nZ
    end  #for iA in 1:nA


    #STOCHASTIC GRID SCHEME....
    println(" Stochastic grid starts ;).... ")
    tol = 1.0e-6
    kiter = 0
    maxdif = 10.0
    N0 = 3
    #dictVF = Dict() #unsure if will be needed....
    global tVF = Vinit
    global itp = interpolate((vGridKComplete, vGridZ, vGridA), Vinit, Gridded(Linear()))
    global tPolicyFnIndx = zeros(Int64,nKComplete,nZ,nA)
    global tPolicyFn = repeat(vGridKComplete,1,nZ,nA)
    global policyfnitp = interpolate((vGridKComplete, vGridZ, vGridA), tPolicyFn, Gridded(Linear()))
    VPrevious = zeros(nKComplete,nZ,nA)
    global vMaxDiff = Float64[]

    while maxdif > tol
        if kiter == 0
            itp = interpolate((vGridKComplete, vGridZ, vGridA), Vinit, Gridded(Linear()))
            policyfnitp = interpolate((vGridKComplete, vGridZ, vGridA), tPolicyFn, Gridded(Linear()))
        end

        # Nkiter = Int( round(2^(1*kiter)*N0 /15) )
        # if Nkiter <3
        #      Nkiter = 3
        # end

        Nkiter = Int( round(2*(1*kiter) /15) )
        if Nkiter <3
            Nkiter = 3
        elseif Nkiter > 43
            Nkiter = Int( round(8*(1*kiter) /15) )
        elseif Nkiter > 20
            Nkiter = Int( round(4*(1*kiter) /15) )
        end

        kiter = kiter+1

        #gets a sample of kapitals...
        Spk = Random.Sampler(MersenneTwister(kiter), 1:(nKComplete))
        rng = MersenneTwister(kiter)
        SkI = rand(rng, Spk,Nkiter,1)  #gets Nkiter-size sample from the possible points
        SkI = sort!(SkI,dims=1) #order it
        SkI[1] = 1
        SkI[end] = nKComplete
        vGridK = vGridKComplete[SkI]

        #       some could be repeated
        #Find the probability matrix P([iZ,iA],ikp,->[ik',iz',ia'])
        PP = zeros(nZ,nA,nKComplete,nZ,nA)
        for iA in 1:nA
            for iZ in 1:nZ
                for ii in SkI

                    suma = 0
                    for zii in 1:nZ
                        for Aii in 1:nA
                            PP[iZ,iA,ii,zii,Aii] = mTranstnZ[iZ, zii] * mTranstnA[iA, Aii]
                            suma = suma + PP[iZ,iA,ii,zii,Aii]
                        end
                    end
                    if suma != 0
                        PP[iZ,iA,ii,:,:] = PP[iZ,iA,ii,:,:]./suma
                    end

                end
            end
        end #loop for PP probabilities
        VPrevious = itp(vGridK[:], vGridZ,vGridA)

        #println("sort of VFI within stochastic starts")
        #update, find max, will require ml1 and ml2
        # CHECK INDICES ARE OK, OW IT WILL DIE....... :O
        #      VFhat, tPolicyFnIndx = gamma(VPrevious)
        VFhat = zeros(Nkiter,nZ,nA)
        tPolicyFnIndx = zeros(Int64,Nkiter,nZ,nA)
        vProvisional = zeros(Nkiter)
        for iA in 1:nA
            for iZ in 1:nZ
                for iK in 1:Nkiter#SkI
                    iKK = SkI[iK]
                    for iKNext in 1:Nkiter#SkI
                        ikkNEXT = SkI[iKNext]
                        #get expected value:
                        expected = 0.0
                        for iAnext in 1:nA
                            for iZnext in 1:nZ
                                expected = expected + PP[ iZ,iA,ikkNEXT,iZnext,iAnext  ] * VPrevious[iKNext,iZnext,iAnext]
                            end
                        end #for iAnext

                        eZ = exp(vGridZ[iZ])
                        L1solved = mL1[iKK,ikkNEXT,iZ,iA]#vGridK[iK]^(1-α)
                        L2solved = mL2[iKK,ikkNEXT,iZ,iA]#vGridK[iK]^(1-α)
                        vProvisional[iKNext] = val_provisional(vGridKComplete[iKK], vGridKComplete[ikkNEXT], vGridZ[iZ],vGridA[iA],L1solved, L2solved, expected, economy)
                    end #for iKNext
                    VFhat[iK, iZ, iA], tPolicyFnIndx[iK, iZ, iA] = findmax(vProvisional)
                end #for iK
            end # for iZ
        end  #end of for iA and fors...

        #maxDiff = norm(VFhat-VPrevious)
        maxdif = maximum(abs.(VFhat-VPrevious))
        vMaxDiff = push!(vMaxDiff,maxdif)

        itp = interpolate((vGridK[:], vGridZ, vGridA), VFhat, Gridded(Linear()))
        tPolicyFn = vGridK[tPolicyFnIndx]
        policyfnitp = interpolate((vGridK[:], vGridZ, vGridA), tPolicyFn, Gridded(Linear()))
        #println("maxdif = ", maxDiff)
        if(mod(kiter,10)==0 || kiter == 1)
            println(" kiter = ", kiter, " Nkiter = ",Nkiter, " Sup Diff = ", maxdif)
        end #if
    end #while maxdif>tol



    tVF = itp(vGridKComplete, vGridZ, vGridA)
    tPolicyFn = policyfnitp(vGridKComplete, vGridZ, vGridA)
    return tVF, tPolicyFn, vGridKComplete, vMaxDiff, mL1, mL2
end

end
