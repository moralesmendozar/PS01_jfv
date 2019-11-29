module test01
export saveTest
using JLD2

function saveTest(doload,dosave,namsave)
  variableTest01 = 0.5
  if doload ==1
    println("loading...")
    #@load namsave  #"test01Data.jld"
    temp = variableTest01 + 1
    println("variableTest01 +1 = ", temp)
  end
  if dosave ==1
    println("saving...")
    @save namsave
    println("variableTest01 = ", variableTest01)
  end
end

end
