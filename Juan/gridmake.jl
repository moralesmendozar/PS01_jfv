module exercise7aux

export gridmake

    """
    GRID CONSTRUCTOR
    Grid with common values of x independently of its size

    Useful for memoization
    """
    function gridmake(x_min::Real, x_max::Real, n::Integer)

        known = Dict(0=>2, 1=>3)

        function mid_points(n)
            if n ∈ keys(known)
                return known[n]
            end
            res = mid_points(n-1) + 2^(n-1)
            known[n] = res
            res
        end
        xgrid = collect(range(x_min, x_max, length = mid_points(n)))

        return xgrid
    end

end
