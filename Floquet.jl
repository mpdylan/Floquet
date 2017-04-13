module Floquet

import ODE.ode23
import ForwardDiff.derivative
importall ApproxFun
import Roots
import Polynomials: Poly, polyval, roots
importall Gadfly
include("OPUC.jl")

const τ = 2π

## Routines for the continuous setting

function truncf(f::Function, truncmax::Float64)
  return x -> (if x < truncmax f(x) else 0.0 end)
end

# The Floquet discriminant for a Schroedinger operator on \R
function disccts(f, λ, P)
  function rhs(t, y)
    return [y[2], f(t)*y[1] - λ*y[1]]
  end
  t1, y1 = ode23(rhs, [1.,0.], [0., P])
  t2, y2 = ode23(rhs, [0.,1.], [0., P])
  0.5*(y1[end][1] + y2[end][2])
end

function perspec(f::Function, P, λmin, λmax)
  Δ(λ) = disccts(f, λ, P)
  
end

function aperspec(f::Function, λmin, λmax)

end

function plotspike(f::Function, P, λmax)
  disc = Fun(λ -> disccts(f, λ, P), Interval(0, λmax))
  heights = [acos(x) - eps()*im for x in roots(disc')]
  locs = [nπ/P for n = 1:length(heights)]
  p = plot(x=locs, xend=locs, y=-1*heights, yend=heights, Geom.segment)
  p
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

function recOPRL(a::Array, b::Array, x::Number)
  p = x
  q = 1.
  tempq = 0.
  for i = 1:length(b)-1
    tempq = p
    p = (1./b[i+1]) * ( (x - a[i])*p - b[i]*q )
    q = tempq
  end
  p
end

function recOPRL2(a::Array, b::Array, x::Number)
  p = 1./b[1]
  q = 0.
  tempq = 0.
  for i = 1:length(b)-1
    tempq = p
    p = (1./b[i+1]) * ( (x - a[i])*p - b[i]*q )
    q = tempq
  end
  p
end

recdiscRL(a::Array, b::Array, x::Number) = recOPRL(a, b, x) - recOPRL2(a[1:end-1], b[1:end-1], x)

function plotspike(a::Array{Float64}, b::Array{Float64})
  disc(x) = recdiscRL(a, b, x)
  ddisc(x) = derivative(disc, x)
  heights = [imag(acos(0.5*disc(x) - eps()*im)) for x in Roots.fzeros(ddisc, -10, 10, no_pts=2000)]
  locations = [n/length(heights) for n=1:length(heights)]
  m = maximum(abs(heights))
  p = plot(x=locations, xend=locations, y=heights, yend=-1*heights, Geom.segment, Coord.cartesian(ymin=-1*m, ymax=m))
  p
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
  A /= (1 - abs(α)^2)
  A
end

# A transfer matrix for Szego recursion a la Simon
function tmatUC(α::Function, z::Number, n::Integer)
  #n = length(α)
  A = eye(2)
  for i = 1:n
    A *= rmatUC(α(i), z)
  end
  A
end

tmatUC(α::Array, z::Number) = tmatUC(i -> α[i], z, length(α))

function opucval(α::Function, z::Number, n::Integer, λ::Number, normed=Bool)
  ϕ = [1. + 0.0im, conj(λ)]
  for i = 1:n
    ϕ = [ϕ[1]*z - conj(α(i)) * ϕ[2], ϕ[2] - α(i)*z*ϕ[1]]
    if normed==true
      ϕ /= sqrt(1 - α(i) * conj(α(i)))
    end
  end
  ϕ .* [1., λ]
end

opucval(α::Array, z::Number, n::Integer, λ::Number, normed=Bool) = opucval(i -> α[i], z, n, λ, normed)

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

function plotspike(α::Array)
  discriminant = Fun(x -> real(discUC(α, exp(im*x))), Interval(0, 2π))
  heights = [imag(acos(0.5*discriminant(x) - eps()*im)) for x in roots(discriminant')]
  locations = [2π*n/length(α) for n=1:length(heights)]
  m = maximum(abs(heights))
  p = plot(x=locations, xend=locations, y=heights, yend=-1*heights, Geom.segment, Coord.cartesian(ymin=-1.5m, ymax=1.5m))
  p
end

function transcoeff(α::Array, p, z)
  ϕ = OPUC.opuc1(α, p)[1][p] / OPUC.opuc1(α, p)[3][p]
  ψ = OPUC.opuc2(α, p)[1][p] / OPUC.opuc2(α, p)[3][p]
  return (z^(-p/2) * (polyval(ϕ, z) + polyval(ψ, z)))
end

function trco(α::Array, z::Number, p::Integer)
  T = tmatUC(α[1:p],z)
  ϕ = (T * [1., 1.])[1]
  ψ = -(T * [1., -1.])[1]
  return z^(-length(α)/2) * (ϕ + ψ)
end

function plotspike_env(α::Array, p)
  discriminant = Fun(x -> real(discUC(α, exp(im*x))), Interval(0, 2π))
  heights = [imag(acos(0.5*discriminant(x) - eps()*im)) for x in roots(discriminant')]
  locations = [2π*n/length(α) for n=1:length(heights)]
#  s(z) = transcoeff(α, p)
  q = plot(layer(x=locations, xend=locations, y=heights, yend=-1*heights, Geom.segment), layer(x -> acosh(0.5abs(transcoeff(α, p, exp(-im*x)))), 0, 2π))
end

function dirichlet(α::Array, z::Number)
  T = tmatUC(α,z)
  ϕ = T * [1., 1.]
  return ϕ[2] - ϕ[1]
end

function dirichletdata(α)
  roots(Fun(x -> dirichlet(α, exp(im*x)), Interval(0, τ)))  
end

function plotspike_dir(α::Array)
  discriminant = Fun(x -> real(discUC(α, exp(im*x))), Interval(0, 2π))
  ϕ, ψ, norms = OPUC.opuc1(α, length(α))
  dirdat = sort(angle(Roots.roots(ψ[end] - ϕ[end]))) + π
  heights = [imag(acos(0.5*real(discUC(α, exp(im*x))) - eps()*im)) for x in roots(discriminant')]

  dirheights = [imag(acos(0.5*real(discUC(α, exp(im*x))) - eps()*im)) for x in dirdat]
#  test = map(>, dirheights, heights)
#  if sum(test) > 0
#    dirheights = vcat(dirheights[end], dirheights[1:end-1])
#  end
  locations = [2π*n/length(α) for n=1:length(heights)]
  m = maximum(abs(heights))
  q = plot(layer(x=locations, xend=locations, y=heights, yend=-1*heights, Geom.segment), layer(x=locations[1:length(dirheights)], y=dirheights, Geom.point))
end


function pspec(α::Array{Complex{Float64}})
  disc = Fun(x -> real(discUC(α, exp(im*x))), Interval(0,2π))
  return roots(disc - 2)
end

function apspec(α::Array{Complex{Float64}})
  disc = Fun(x -> real(discUC(α, exp(im*x))), Interval(0,2π))
  return roots(disc + 2)
end

function spiketop(α::Array{Complex{Float64}})
  disc = Fun(x -> real(discUC(α, exp(im*x))), Interval(0,2π))
  return roots(disc')
end

#function pspec(f::Function, P::Float64, max::Float64)
#  disc = Fun(λ -> disccts(f, λ, P), 0..max)
#  return roots(disc - 1)
#end

#function apspec(f::Function, P::Float64, max::Float64)
#  disc = Fun(λ -> disccts(f, λ, P), 0..max)
#  return roots(disc + 1)
#end

#function spiketop(f::Function, P::Float64, max::Float64)
#  disc = Fun(λ -> disccts(f, λ, P), 0..max)
#  return roots(disc')
#end



#function spikepeaks(α::Array{Complex{Float64}})
#
#end
#
#function spikepeaks(f::Function, P::Float64, max::Float64)
#
#end


function decomp(α, z, sign)
  A = tmatUC(α, z)
  eigval, evec = eig(A)
  evec = evec[:, sortperm(eigval, by=abs)]
  v = evec[:,sign]
  η = (v[1] + v[2]) / 2
  ζ = (v[1] - v[2]) / 2
  (η, ζ)
end

#=
function vbaker(α::Function, p::Integer, n::Integer, z::Number, sign::Integer)
  coeff(k::Integer) = α[(x % p) + 1]
  Δ = discUC(α, z)
  ϕn = opucval(α, z, n, 1.)
  ϕp = opucval(α, z, p, 1.)
  ψn = opucval(α, z, n, -1.)
  ψp = opucval(α, z, p, -1.)
  βn = Δ - 2ψn + sign*sqrt(Δ^2 - 4)
  ηn = -2ϕn - Δ - sign*sqrt(Δ^2 - 4)
  βp = Δ - 2ψp + sign*sqrt(Δ^2 - 4)
  ηp = -2ϕn - Δ - sign*sqrt(Δ^2 - 4)
end

function baker(α::Array, n::Integer, z::Number, sign::Integer)
  coeff(k::Integer) = α[(k % length(α)) + 1]
  A = tmatUC(α, z)
  eigval, evec = eig(A)
  if abs(trace(z^(-length(α) / 2) * A)) > 2
    evec = evec[:, sortperm(eigval, by=abs)]
  end
  denom = evec[:, sign]
  num = evec[:, sign]
  for i = 1:n
    num = rmatUC(coeff(i), z) * num
  end
  return (num[1] - num[2]) / (denom[1] - denom[2])
end

function adjbaker(α::Array, n::Integer, z::Number, sign::Integer)
  coeff(k::Integer) = α[(k % length(α)) + 1]
  A = tmatUC(α, z)
  eigval, evec = eig(A)
  if abs(trace(z^(-length(α) / 2) * A)) > 2
    p = sortperm(eigval, by=abs)
    evec = evec[:, p]
    permute!(eigval, p)
  end
  denom = evec[:, sign]
  num = evec[:, sign]
  for i = 1:n
    num = rmatUC(coeff(i), z) * num
  end
  return eigval[sign]^(-n / length(α)) * (num[1] - num[2]) / (denom[1] - denom[2])
end
=#

function normbaker(α::Array, n::Integer, z::Number, sign::Integer)
  z^(-length(α) / 2) * adjbaker(α, n, z, sign)
end

function baker2(α::Array, n::Integer, z::Number, sign::Integer)
  coeff(k::Integer) = α[(k % length(α)) + 1]
  A = tmatUC(α, z)
  eigval, evec = eig(A)
  if abs(trace(z^(-length(α) / 2) * A)) > 2
    evec = evec[:, sortperm(eigval, by=abs)]
  end
  denom = evec[:, sign]
  num = evec[:, sign]
  for i = 1:n
    num = rmatUC(coeff(i), z) * num
  end
  return (num[1] + num[2]) / (denom[1] + denom[2])
end


end
