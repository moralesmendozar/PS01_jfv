module ex7aux

export gridNum

    """
    GRID CONSTRUCTOR
    Grid with common values of x independently of its size

    Useful for memoization
    """
    function gridNum(x_min::Real, x_max::Real, n::Integer)

        #known = Dict(0=>2, 1=>3)

        #function mid_points(n)
    #        if n ∈ keys(known)
    #            return known[n]
    #        end
    #        res = mid_points(n-1) + 2^(n-1)
    #        known[n] = res
    #        res
    #    end
        #xgrid = collect(range(x_min, x_max, length = mid_points(n)))

        function almostFib(n)
                if n == 1
                    return 3
                else
                    res = almostFib(n-1) + 2^(n-1)
                    return res
                end
        end

        if n<1
            println("Number n has to be positive!!")
        end

        xgrid = collect(range(x_min, x_max, length = almostFib(n)))

        return xgrid
    end

end
