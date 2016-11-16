module Floquet

import ODE.ode23
importall Polynomials
#importall OPUC

function spike(disc::Function, domain, p::Number)
  #δ = domain[2] - domain[1]
  y = map(disc, domain) + 0.*im
  #dy = diff(y) / δ
  dy = diff(y)
  spre = zeros(dy)
  spim = zeros(dy)
  spsq = sqrt(1 - y .^ 2)
  for i = 1:(length(dy)-1)
    spre[i+1] = spre[i] + abs(real(dy[i] / spsq[i]))
    spim[i+1] = spim[i] + imag(dy[i] / spsq[i])
  end
  return (1/p)*spre, spim
end

## Routines for the continuous setting

# The Floquet discriminant for a Schroedinger operator on \R
function disccts(f, λ, P)
  function rhs(t, y)
    return [y[2], f(t)*y[1] - λ*y[1]]
  end
  t1, y1 = ode23(rhs, [1.,0.], [0., P])
  t2, y2 = ode23(rhs, [0.,1.], [0., P])
  0.5*(y1[end][1] + y2[end][2])
end

## Routines for Jacobi matrices (OPRL)

function permat(pot::Array, h)
  n = 2*length(pot)
  C = zeros(typeof(h), n, n)
  for i = 1:n-1 
    C[i, i+1] = 1 
  end
  C += transpose(C)
  for i = 1:div(n,2)
    C[i, i] = pot[i]
    C[i + div(n,2), i + div(n,2)] = pot[i]
  end
  C[1, n] = h
  C[n, 1] = 1/h
  return C
end

# A single recursion matrix for OPRL/Jacobi
function rmatRL(a::Number, x::Number)
  A = Array(Float64, 2, 2)
  A[1,1] = x - a
  A[1,2] = -1
  A[2,1] = 1
  A[2,2] = 0
  A
end

# Transfer matrix for OPRL/Jacobi
function tmatRL(a::Array, x::Number)
  A = eye(2)
  for i = 1:length(a)
    A *= rmatRL(a[i], x)
  end
  A
end

# Floquet discriminant for OPRL/Jacobi 
function discRL(pot::Array, x::Number)
  n = length(pot)
  A = SymTridiagonal(pot, ones(typeof(pot[1]),n-1)) - SymTridiagonal(fill(x, n), zeros(n-1))
  B = SymTridiagonal(pot[2:n-1], ones(typeof(pot[1]), n-3)) - SymTridiagonal(fill(x, n-2), zeros(n-3))
  return (det(A) - det(B))
end

# Incomplete, fix this
function rmatGC(α::Number, z::Number, j::Int64)
  A = Array(typeof(α), 2, 2)
  A[1,1] = z
  A[1,2] = -1 * conj(α) * z^j
  A[2,1] = -α * z^(-j)
  A[2,2] = z^-1
  A /=  #Need to fill in the correct norming constant here!
  A
end

## Routines for Szego recursion (OPUC)

# A CMV matrix, periodized

function periodCMV(α::Array, n::Integer, β::Number)
  M = [-α[n]]
  ρ = sqrt(1 - abs(α[1])^2)
  L = reshape([conj(α[1]), ρ, ρ, -α[1]], 2, 2)
  for i = 2:n-1
    ρ = sqrt(1 - abs(α[i])^2)
    θ = reshape([conj(α[i]), ρ, ρ, -α[i]], 2, 2)
    if iseven(i-1)
      M = cat([1,2], M, θ)
    else
      L = cat([1,2], L, θ)
    end
  end
  M = cat([1,2], M, [conj(α[n])])
  M[1,n] = (1/β) * sqrt(1 - abs(α[1])^2)
  M[n,1] = β * sqrt(1 - abs(α[n])^2)
  L, M
end

# A recurrence matrix
function rmatUC(α::Number, z::Number)
  A = Array(typeof(α), 2, 2)
  A[1,1] = z
  A[1,2] = -1 * conj(α)
  A[2,1] = -α * z
  A[2,2] = 1
  A /= sqrt(1 - α*conj(α))
  A
end

# A transfer matrix for Szego recursion a la Simon
function tmatUC(α::Array, z::Number)
  n = length(α)
  A = eye(2)
  for i = 1:n
    A *= rmatUC(α[i], z)
  end
  A
end

# The transfer matrix for the Geronimo-Case version of OPUC recursion
function tmatGC(α::Array, z::Number)
  A = tmatUC(α, z)
  A *= diagm([1, z^(-n)])
end

function discUC(α::Array, z::Number)
  n = length(α)
  if isodd(n)
    error("For annoying technical reasons we assume even period. Try a vector of even length")
  end
  return z^(-n/2) * trace(tmatUC(α, z))
end

function trco(α::Array, z::Number)
  T = tmatUC(α,z)
  ϕ = (T * [1., 1.])[1]
  ψ = -(T * [1., -1.])[1]
  return z^(-length(α)/2) * (ϕ + ψ)
end

function dirichlet(α::Array, z::Number)
  T = tmatUC(α,z)
  ϕ = T * [1., 1.]
  return ϕ[2] - ϕ[1]
end

end
