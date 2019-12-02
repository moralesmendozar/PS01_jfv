#interpolate learning...
using Interpolations, NLsolve

# Now learning to do interpolation for 2D or so...
println("Interpolation learning stuff...")

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
#itp = interpolate(gridK,VF[:,1,1],Gridded(Linear()))
#itp(gridK[[1,4,10]])
itp = interpolate((gridK,vZ,vA),VF,Gridded(Linear()))


gridK2 = gridK[[1,3,5,5,7,10,15,20]]
itp[gridK2, vZ,vA]

#pick the random from 1:nK, 1:nA,
nKsample = 5
nAsample = 2
nZsample = 4
vkRandChoices = sort(shuffle(Vector(1:nK))[1:nAsample])
vARandChoices = sort(shuffle(Vector(1:nA))[1:nAsample])
vZRandChoices = sort(shuffle(Vector(1:nZ))[1:nZsample])
