module exercise2

export ss_solve

using Parameters, NLsolve

function ss_solve(economy, initialGuess::Vector)

    @unpack α, β, θ, δ = economy

    function foc!(F,x)

        c1 = x[3]^α * x[1]^(1-α) - δ * x[3]

        # FOC w.r.t labor 1
        F[1] = θ * (1-α)*x[3]^α*x[1]^(-α) * (c1)^(θ-1)*x[2]^(1-θ) - x[1] - x[2]

        # FOC w.r.t labor 2
        F[2] = (1 - θ)* (c1)^θ * x[2]^(-θ) - x[1] - x[2]

        # FOC w.r.t. capital
        F[3] =  β* α * (x[1])^(1-α) * (x[3])^(α-1) +  β*(1 - δ) - 1

    end

    function soc!(J,x)

        c1 = x[3]^α * x[1]^(1-α) - δ * x[3]

        # Derivatives of the 1st equation
        aux1 = θ * (1-α) * x[3]^(α) * x[1]^(-α)

        J[1,1] = aux1 * ( -α *(c1)^(θ-1)*x[2]^(1-θ)* x[1]^(-1) + (θ-1)* c1^(θ-2)* x[2]^(1-θ) * (1-α)*x[3]^(α) * x[1]^(-α)   ) - 1

        J[1,2] = θ*(1-θ)* c1^(θ-1) * x[2]^(-θ)* (1-α) * x[3]^α * x[1]^(-α) - 1

        J[1,3] = aux1 * ( α * (c1)^(θ - 1) *x[2]^(1-θ) * x[3]^(-1)    +    (θ-1) *(c1)^(θ-2)*x[2]^(1-θ)* (α*x[3]^(α-1)*x[1]^(1-α) - δ)  )

        # Derivatives of the 2nd equation
        J[2,1] = θ * (1-θ )* c1^(θ - 1) * x[2]^(-θ) * (1-α) * x[3]^(α) * x[1]^(-α)  - 1

        J[2,2] = -θ * (1-θ ) * c1^(θ) * x[2]^(-1-θ) - 1

        J[2,3] = θ * (1-θ )* c1^(θ - 1) * x[2]^(-θ) * ( α * x[3]^(α-1) * x[1]^(1-α)  - δ)

        # Derivatives of the 3rd equation
        J[3,1] = α *β *(1 - α) * x[3]^(α-1) * x[1]^(-α)

        J[3,2] = 0

        J[3,3] = α *β *(α - 1) * x[3]^(α-2) * x[1]^(1-α)

    end

    results = nlsolve(foc!,soc!, initialGuess)

    return results
end


end
