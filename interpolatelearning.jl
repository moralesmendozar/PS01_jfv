#interpolate learning...
using Interpolations, NLsolve
#Learning to do random stuff first:
rng = MersenneTwister(2);
xs = shuffle(rng, Vector(1:10))

Xs = xs[1:6]

println("Xs = ", Xs')
println("done with xs, now do ns")

rng = MersenneTwister(123)
sp = Random.Sampler(rng, 1:10) # or Random.Sampler(MersenneTwister, 1:20)
Skk = rand(rng, sp,20,1)
println("Skk = ", Skk)
sort!(Skk,dims=1)
println("sorted! :) ", Skk)
for ii in 1:20
    n = rand(rng, sp) # similar to n = rand(rng, 1:20)
    #println("random n = ", n)
    # use n
end

vectorRand = rand(rng, sp,1,20)
println("vectorRand = ", vectorRand)
println("done with Random, now going into Interpolation stuff...")

# Now learning to do interpolation for 2D or so...

gridK = range(1,stop=10,length=20)
nK = length(gridK)
vZ    = [-0.0673, -0.0336,  0,      0.0336, 0.0673]
vA    = [  0.9,  1,    1.1]
nZ = 5
nA = 3
VF = repeat(gridK,1,nZ,nA)

for iA in 1:nA
    for iZ in 1:nZ
        for iK in 1:nK
            VF[iK,iZ,iA] = exp(vZ[iZ])*exp(vA[iA])*VF[iK,iZ,iA]
        end
    end
end

#INTERPOLATION LEARNING:
interpolate(gridK,VF[:,1,1],Gridded(Linear()))



#pick the random from 1:nK, 1:nA,
nKsample = 5
nAsample = 2
nZsample = 4
vkRandChoices = sort(shuffle(Vector(1:nK))[1:nAsample])
vARandChoices = sort(shuffle(Vector(1:nA))[1:nAsample])
vZRandChoices = sort(shuffle(Vector(1:nZ))[1:nZsample])

println("Quick check on Nk = 2^(2*kiter)*N0")
N0 = 10
println("N0 = ", N0)
for kiter in 0:10
    Nk = 2^(2*kiter)*N0
    println("2^(2*", kiter,")*N0 = ", Nk)
end
