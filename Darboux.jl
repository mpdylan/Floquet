module Darboux

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

dl(n) = eye(n) + -1*rshift(n)
dr(n) = -1*eye(n) + 1*lshift(n)

# Darboux transformation functions

function factor_reverse(λ, u)
    μl = ( dl(n) * u ) ./ u
    μr = ( dr(n) * u ) ./ u
    (dr(n) + diagm(ul)) * (dl(n) + lshift(n)*diagm(ur)) + λ*eye(n)
end
