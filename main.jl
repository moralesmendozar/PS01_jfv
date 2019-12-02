# UPENN
# 714, computational econ
# profr. JssFdzVllvrd
# PS01
# Rodrigo A Morales M
# November, 2019
# ------------------------------------------------------------------------------
# 0. init, packages, and modules
println("Retrieving packages...")
using Parameters, Plots, JLD2, JLD
println("Packages gotten alright, getting modules...")
cd("/home/moralesmendozar/Dropbox/03_UPenn/classes/2019_fall/00_714/partJesus/PS01_jfv")
include("steadyState.jl")
using .ss_deterministic
include("a_fixed_grid.jl")
using .ex3
include("a_fixed_grid_optimized.jl")
using .ex3b
include("b_accelerator.jl")
using .ex03b
include("c_multigrid.jl")
using .ex03c
include("c_multigrid_enhanced.jl")
using .ex3c
include("d_stochastic.jl")
using .ex3d
include("egm_general.jl")
using .ex4
println(" modules loaded ok")
# ------------------------------------------------------------------------------
# 1. Set Parameters (building structure)
econ_params = @with_kw (
                α = 0.33,   # kaptl share
                β = 0.96,   # time discount
                δ = 0.1,    # deprectn
                θ = 0.5,    # elast subs tween goods
                maxiter = 1000,
                # Stochasticity
                vGridZ    = [-0.0673, -0.0336,  0,      0.0336, 0.0673],
                mTranstnZ = [ 0.9727  0.0273    0       0       0;
                              0.0041  0.9806    0.0153  0       0;
                              0       0.0082    0.9836  0.0082  0;
                              0       0         0.0153  0.9806  0.0041;
                              0       0         0       0.0273  0.9727],
                # production func. 2
                vGridA    = [  0.9,  1,    1.1],
                mTranstnA = [  0.9   0.1   0;
                               0.05  0.9   0.05;
                               0     0.1   0.9],
                mTranstnZA = kron(mTranstnZ,mTranstnA)
        )
# once structure is built, call it to use in VFIs (fixed grid, etc...)
econparams = econ_params()  #call function and get the relevant parameters 4 ss
dosave = 0
doload = 0
@unpack α, β, δ, θ, vGridZ, vGridA  = econparams
# ------------------------------------------------------------------------------
# 2. Steady State, solve for:
xinit =  [0.5, 0.5, 1] #[0.2, 0.2, 0.8] #[0.5, 0.5, 1] #initial value to find ss
sssol =  ss_deterministic.ss_solver(xinit, α, β, δ, θ )
l1ss  =  sssol.zero[1]
l2ss  =  sssol.zero[2]
kss   =  sssol.zero[3]

# the rest are determined by k,l1,l2
lss   =  l1ss + l2ss
c1ss  =  kss^(α) * l1ss^(1-α) - δ * kss
c2ss  =  l2ss

println(" Steady state ")
println("K_ss = ", kss, "  L_ss = ", lss, "  C1_ss = ", c1ss, " C2_ss = ", c2ss)
println("---------------------------------------------------------------------")
println(" ")

SteadyState_Variables = @with_kw (
                        kss  = kss,
                        l1ss =  l1ss,
                        l2ss =  l2ss,
                        lss  =  lss,
                        c1ss =  c1ss,
                        c2ss =  c2ss,
                        utilitySS = c1ss^(0.5) * c2ss^(0.5) - 0.5*(lss)^2
)
SSVarbls = SteadyState_Variables()
# ------------------------------------------------------------------------------
# 03.  PS 01.ex3) Fixed grid:
if doload == 1
        #@load "tmpfile.jld"
        #mVinit = mVF
else
        @unpack kss, utilitySS= SSVarbls
        nk = 250
        nZ = length(vGridZ)
        nA = length(vGridA)
        vGridK = collect(range(0.7 * kss, 1.3 * kss, length = nk))
        mVinit = repeat(vGridK,1, nZ,nA)
        #mVinit = fill(utilitySS, nk, nZ, nA)
end

## FIXED GRID
#println(" calling a_fixed_grid_optimized... ")
#@time mVF, mPolicyFn, vGridK = ex3.a_fixed_grid(econparams, SSVarbls,mVinit,nk)
#@time mVF, mPolicyFn, vGridK, vMaxDifference = ex3b.a_fixed_grid_optimized(econparams, SSVarbls,mVinit,nk,0,0,"OptimFunData.jld")

## ACCELERATOR
#println(" calling accelerator... ")
#@time mVF, mPolicyFn, vGridK, vMaxDifference = ex03b.b_accelerator(econparams, SSVarbls,mVinit,nk)

##  MultiGrid:
#println(" calling multigrid... ")
##   f(1) = 3 f(2) = 5 f(3) = 9 f(4) = 17 f(5) = 33 f(6) = 65 f(7) = 129 f(8) = 257
##   f(9) = 513 f(10) = 1025 f(11) = 2049 f(12) = 4097 f(13) = 8193 f(14) = 16385
##   f(15) = 32769 f(16) = 65537 f(17) = 131073 f(18) = 262145
##   f(19) = 524289   f(20) = 1048577
nMidPoints = [4 5] ## [4 5 7] easy   #  real stuff: [5 7 9] [5, 7, 9, 12]
#@time mVF, mPolicyFn, vGridK, vMaxDifference = ex03c.c_multigrid(econparams, SSVarbls, nMidPoints)
#@time mVF, mPolicyFn, vGridK, vMaxDifference = ex3c.c_multigrid_enhanced(econparams, SSVarbls, nMidPoints)

## ENDOGENOUS
println(" calling endogenous Grid... ")
@time tValueFunctionTilde = ex4.egm_general(econparams, SSVarbls, nk)
println("using endogenous as input for accelerator:")
@time mVF, mPolicyFn, vGridK, vMaxDifference = ex03b.b_accelerator(econparams, SSVarbls,tValueFunctionTilde,nk)

#STOCHASTIC METHOD IS A BIT TRICKIER...
# Grid Capital
nKComplete = 1000
vGridKComplete = collect(range(0.85 * kss,1.25 * kss, length = nKComplete))
###mVinitStoch = repeat(vGridKComplete,1, nZ,nA)
mVinitStoch = fill(utilitySS, nKComplete, nZ, nA)
#@time mVF, mPolicyFn, vGridK, vMaxDifference, mL1, mL2 = ex3d.d_stochastic(econparams, SSVarbls,vGridKComplete, nKComplete,mVinitStoch)
# ------------------------------------------------------------------------------
# 04. Graphs 03.03 Value Functions:
pValueFunction = plot(vGridK, mVF[:,1,1],title="Value Function", label = "z_1, A_1", xlabel = "Capital",legend=:topleft)
plot!(vGridK, mVF[:,end,1], label = "z_5, A_1")
plot!(vGridK, mVF[:,1,end], label = "z_1, A_5")
plot!(vGridK, mVF[:,end,end], label = "z_5, A_3")
#Save plots...
savefig("Plots/003_ValueFunction_endogenous50_20191201.png")
# Graphs Policy Functions:
pPolicyFunction =  plot(vGridK,vGridK,title="Policy Function", color=:black,linestyle=:dash)
plot!(vGridK, mPolicyFn[:,1,1], label = "z_1, A_1", xlabel = "Capital", color=:blue)
plot!(vGridK, mPolicyFn[:,end,1],label = "z_5, A_1", color=:blue, linestyle=:dash)
plot!(vGridK, mPolicyFn[:,1,end], label = "z_1, A_5", color=:red)
plot!(vGridK, mPolicyFn[:,end,end], color=:red, linestyle=:dash, label = "z_5, A_3")
# savefig
savefig("Plots/003_PolicyFunction_endogenous50_20191201.png")
#plot(pValueFunction)
#plot(pPolicyFunction)



# ------------------------------------------------------------------------------
#       a_fixed_grid does:
#Iteration = 360 Sup Diff = 9.836069397986712e-7
#23345.156622 seconds (15.06 G allocations: 260.431 GiB, 0.21% gc time)
#       a_fixed_grid_optimized does:
#Iteration = 360 Sup Diff = 9.836069397986712e-7
#329.447920 seconds (9.70 G allocations: 144.820 GiB, 4.07% gc time)
#       b_accelerator does:
#Iteration = 360 Sup Diff = 9.825041917304719e-7
# 28.158447 seconds (1.04 G allocations: 15.474 GiB, 7.39% gc time)

#       multigrid does
#               for: nMidPoints = [4 5 7]
# Iteration = 70 Sup Diff = 9.803856641972994e-7
# 82.542657 seconds (1.98 G allocations: 31.296 GiB, 7.88% gc time)
#               for: nMidPoints = [5 7 9]
#Iteration = 30 Sup Diff = 1.2622364240652081e-6
# 923.525038 seconds (16.66 G allocations: 263.848 GiB, 9.52% gc time)

#       multigrid_enhanced does
#               for nMidPoints = [4 5]
#Iteration = 50 Sup Diff = 9.855368940455396e-7
#  7.445573 seconds (183.23 M allocations: 3.274 GiB, 7.74% gc time)
#               for: nMidPoints = [4 5 7]
#Iteration = 70 Sup Diff = 9.803856641972994e-7
# 62.945592 seconds (1.98 G allocations: 31.257 GiB, 15.73% gc time)
#               for: nMidPoints = [5 7 9]
#Iteration = 30 Sup Diff = 1.2622364240652081e-6
#507.859205 seconds (16.65 G allocations: 263.614 GiB, 15.72% gc time)
#               for: nMidPoints = [5 7 9 12]

#       stochastic does:
#kiter = 320 Nkiter = 85 Sup Diff = 2.2284117498772016e-6
# kiter = 330 Nkiter = 175 Sup Diff = 1.5725093996388217e-6
# 43.567893 seconds (1.51 G allocations: 23.395 GiB, 10.46% gc time)

#       Endogenous with accelerator does:
#Iteration = 470 Sup Diff = 1.3941744162701624e-7
# 1.323211 seconds (1.53 M allocations: 1.267 GiB, 14.75% gc time)
#using endogenous as input for accelerator:
#and accelerator Iteration = 390 Sup Diff = 1.0642538354482332e-6
# 30.165119 seconds (1.15 G allocations: 17.108 GiB, 7.69% gc time)


if dosave == 1
        #@save "Data/endogenous_20191201.jld"
end
