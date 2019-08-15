# Commutation methods for spectral deformations of Jacobi operators

module Commutation

# Basic matrices

Δ(n) = SymTridiagonal(fill(-2.0, n), ones(n-1))
Δ0(n) = SymTridiagonal(zeros(n), ones(n-1))

function lshift(n)
    A = circshift(eye(n), [0,1])
    A[end,1] -= 1.0
    A
end

function rshift(n)
    A = circshift(eye(n), [0,-1])
    A[1,end] -= 1.0
    A
end

dl(n::Integer) = eye(n) + -1*rshift(n)
dr(n::Integer) = -1*eye(n) + 1*lshift(n)

# Darboux transformation functions, functional setting
dl(f::Function) = (n -> f(n) - f(n - 1))
dr(f::Function) = (n -> f(n + 1) - f(n))

μl(f::Function) = (n -> (f(n) - f(n - 1)) / f(n))
μr(f::Function) = (n -> (f(n + 1) - f(n)) / f(n))

ushift(f::Function) = (n -> dl(μr(f))(n) + dr(μl(f))(n))

# Darboux transformation functions

#function factor_reverse(λ::Number, u::Array)
#    n = length(u)
#    μl = ( dl(n) * u ) ./ u
#    μr = ( dr(n) * u ) ./ u
#    (dr(n) + diagm(μl)) * (dl(n) - rshift(n)*diagm(μr)) + λ*eye(n)
#end

# Compute a Weyl function using zero initial condition at either end
function weyl(J::Array{Float64, 2}, λ::Number, condition::Symbol)
    if condition == :-
        result = Array{Float64}(size(J)[1])
        result[1] = eps()
        result[2] = eps()
        for i = 3:(length(result))
            result[i] = (λ - J[i-1,i-1]) * result[i-1] - result[i-2]
        end
    end
    if condition == :+
        result = Array{Float64}(size(J)[1])
        result[end] = eps()
        result[end-1] = eps()
        for i = length(result)-2:-1:1
            result[i] = (λ - J[i+1,i+1]) * result[i+1] - result[i+2]
        end
    end
    result
end

# Perform single commutation using an "eigenvalue" and a Weyl function
function factor_reverse(u::Array{Float64, 1}, λ::Number)
    n = length(u)
    μl = ( dl(n) * u ) ./ u
    μr = ( dr(n) * u ) ./ u
    (dl(n) - rshift(n)*diagm(μr)) * (dr(n) + lshift(n)*diagm(μl)*rshift(n)) + λ*eye(n)
end

function factor_reverse(J::Array{Float64, 2}, λ::Number, condition::Symbol)
    u = weyl(J, λ, condition)
    factor_reverse(u, λ)
end

# Perform single commutation using the formulas from Teschl

function offdiag(n::Integer, b::Function, u::Function)
    -sqrt( b(n) * b(n + 1) * u(n) * u(n + 2) ) / u(n + 1)
end
function diagentry(n::Integer, a::Function, u::Function, λ::Number)
    λ - a(n) * (u(n) / u(n+1) + u(n+1) / u(n))
end

end
