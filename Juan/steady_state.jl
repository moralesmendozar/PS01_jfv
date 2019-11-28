module exercise2

export ss_solve

using NLsolve

function ss_solve()

    function f!(F,x)
        α = 0.33; β = 0.96; δ = 0.9

        c1 = x[1]^α * x[2]^(1-α) + δ*x[1] - x[1]

        # FOC w.r.t. labor 1
        F[1] = 0.5 * (1-α) * x[1]^α * x[2]^(-α) * (c1)^(-0.5) * x[3]^(0.5) - x[2] - x[3]

        # FOC w.r.t. labor 2
        F[2] = 0.5 * (c1)^(0.5) * x[3]^(-0.5) - x[2] - x[3]

        # FOC w.r.t. capital
        F[3] =  α * β * x[1]^(α-1) * x[2]^(1-α) + β * δ - 1

    end

    function j!(J,x)
        α = 0.33; β = 0.96; δ = 0.9
        c1 = x[1]^α * x[2]^(1-α) + δ*x[1] -x[1]

        # Derivatives of the 1st equation
        aux1 = 0.5 * (1-α) * x[2]^(-α) * x[3]^(0.5)
        J[1,1] = aux1 * (α * x[1]^(α-1)*(c1)^(-0.5) - 0.5 * x[1]^α * (α*x[1]^(α-1)*x[2]^(1-α) + δ - 1) * (c1)^(-1.5) )

        aux2 = 0.5 * (1-α) * x[1]^(α) * x[3]^(0.5)
        J[1,2] = aux2 * ( -α * x[2]^(-α-1)*(c1)^0.5 - x[2]^(-α) * ( 0.5 * (1-α)*x[1]^(α) * x[2]^(-α) * c1^(-1.5) ) ) - 1

        J[1,3] = 0.5 * (1-α) * x[1]^α * x[2]^(-α)* c1^(-0.5) * 0.5 * x[3]^(-0.5) - 1


        # Derivatives of the 2nd equation
        J[2,1] = 0.5 * x[3]^(-0.5) * (0.5 * (α * x[1]^(α-1) * x[2]^(1-α) + δ - 1) * c1^(-0.5))

        J[2,2] = 0.5 * x[3]^(-0.5) * (0.5 * (1-α) * x[1]^(α) * x[2]^(-α) * c1^(-0.5) ) - 1

        J[2,3] = -0.5 * c1^(0.5) * x[2]^(-1.5) - 1


        # Derivatives of the 3rd equation
        J[3,1] = α * (α - 1) * x[1]^(α-2) * x[2]^(1-α)

        J[3,2] = α * x[1]^(α-1) * (1-α) * x[2]^(-α)

        J[3,3] = 0

    end

    results = nlsolve(f!,j!, [0.8, 0.2, 0.2])

    return results
end


end
