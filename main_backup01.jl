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
println("Packages gotten alright.")
cd("/home/moralesmendozar/Dropbox/03_UPenn/classes/2019_fall/00_714/partJesus/PS01_jfv")
include("steadyState.jl")
using .ss_deterministic
println("SteadyState module ok")
include("a_fixed_grid.jl")
using .ex3
println("Fixed Grid module ok")
include("a_fixed_grid_optimized.jl")
using .ex3b
println("Fixed Grid Optim module ok")
include("b_accelerator.jl")
using .ex03b
println("b_accelerator module ok")
include("saveTest.jl")
#using .test01
#println("Save Test module ok")
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
        @load "tmpfile.jld"
        mVinit = mVF
else
        @unpack kss = SSVarbls
        nk = 250
        nZ = length(vGridZ)
        nA = length(vGridA)
        vGridK = collect(range(0.7 * kss, 1.3 * kss, length = nk))
        mVinit = repeat(vGridK,1, nZ,nA)
end
println(" calling fixed grid optimized solution (03)... ")
#@time mVF, mPolicyFn, vGridK = ex3.a_fixed_grid(econparams, SSVarbls,mVinit,nk)
@time mVF, mPolicyFn, vGridK = ex3b.a_fixed_grid_optimized(econparams, SSVarbls,mVinit,nk,0,0,"OptimFunData.jld")
#test01.saveTest(0,1,"OptimFunData.jld")
# ------------------------------------------------------------------------------
# 04. Graphs 03.03 Value Functions:
pValueFunction = plot(vGridK, mVF[:,1,1],title="Value Function", label = "z_1, A_1", xlabel = "Capital",legend=:topleft)
plot!(vGridK, mVF[:,end,1], label = "z_5, A_1")
plot!(vGridK, mVF[:,1,end], label = "z_1, A_5")
plot!(vGridK, mVF[:,end,end], label = "z_5, A_3")
#Save plots...
savefig("Plots/001_ValueFunction_a_fixed_grid_optimized_20191129.png")
# Graphs Policy Functions:
pPolicyFunction =  plot(vGridK,vGridK,title="Policy Function", color=:black,linestyle=:dash)
plot!(vGridK, mPolicyFn[:,1,1], label = "z_1, A_1", xlabel = "Capital", color=:blue)
plot!(vGridK, mPolicyFn[:,end,1],label = "z_5, A_1", color=:blue, linestyle=:dash)
plot!(vGridK, mPolicyFn[:,1,end], label = "z_1, A_5", color=:red)
plot!(vGridK, mPolicyFn[:,end,end], color=:red, linestyle=:dash, label = "z_5, A_3")
#savefig("Plots/001_PolicyFunction_a_fixed_grid_optimized_20191129.png")
#plot(pValueFunction)
#plot(pPolicyFunction)



# ------------------------------------------------------------------------------
#a_fixed_grid_optimized does:
#Iteration = 360 Sup Diff = 9.836069397986712e-7
#329.447920 seconds (9.70 G allocations: 144.820 GiB, 4.07% gc time)


if dosave == 1
        @save "Data/a_fixed_grid_optimized__20191129.jld"
end
