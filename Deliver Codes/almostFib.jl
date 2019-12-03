function almostFib(n)
        if n == 1
            return 3
        else
            res = almostFib(n-1) + 2^(n-1)
            return res
        end
end

for ii in 1:20
    #println("almostFib(",ii,") = ", almostFib(ii))
    println("f(",ii,") = ", almostFib(ii))
end
