#############################################################################
# Types and functions for analysis of orthogonal polynomials
#############################################################################

module OPUC

importall Polynomials

export opuc, opuc1, opuc2

const τ = 2π

const Z = Poly([0.,1. + 0.*im])
const X = Poly([0.,1.])

#The following functions generate orthogonal polynomials using the recurrence relations.
#Recurrence coefficients are passed as functions. If you have a sequence of recurrence coefficients, pass the anonymous function n -> alpha(n).

#sqrt(f::Function) = x -> sqrt(f(x))
#square(f::Function) = x -> f(x)^2

#One recursive step in Szego recursion
function opuc_recurse(α::Number, ϕ::Poly, ψ::Poly)
  return (ϕ*Z - conj(α)*ψ, ψ - α*Z*ϕ)
end

#Generate the first n orthogonal polynomials on the unit circle based on Verblunsky coefficients alpha.
#λ is an initial condition, and may be any value on the unit circle. If set to 1, ϕ is the orthogonal polynomial and ψ is the reversed polynomial.
#If λ = -1, then ψ becomes the reversed second kind polynomial.
#If λ is another point on the circle, each is a linear combination of first and second kind polynomials.
#If λ is not on the unit circle, you fucked up.

function opuc(α::Function, n::Integer, λ::Complex)
  ϕ = Array(Poly{Complex{Float64}}, n)
  ψ = Array(Poly{Complex{Float64}}, n)
  norms = Array(Float64, n)
  norms[1] = 1
  (ϕ[1], ψ[1]) = opuc_recurse(α(1), one(Poly{Complex{Float64}}), Poly([conj(λ)]))

#  ϕ[1] = one(Poly{Complex{Float64}})
#  ψ[1] = Poly([conj(λ)])
  for j = 2:n
    (ϕ[j], ψ[j]) = opuc_recurse(α(j),ϕ[j-1], ψ[j-1])
  end
  for j = 2:n
    norms[j] = norms[j-1]*(1 - abs(α(j))^2)
  end
  return (ϕ, λ * ψ, norms)
end

opuc(α::Array, n::Integer, λ::Complex) = opuc(i -> α[i], n, λ)

opuc1(α, n) = opuc(α, n, 1. + 0.0im)
opuc2(α, n) = opuc(α, n, -1. + 0.0im)

#Functions for approximating measures on the real line and unit circle
#On the real line we use Edelman's nearly-Toeplitz approximation to obtain a continuous measure
#There is an alternative which returns a discrete measure by using the quadrature formula which obtains the weights by evaluating second kind polynomials.
#On the unit circle we use Bernstein-Szego approximations.
#There is no discrete version of OPUC measures--OPUC zeros lie in the disk, so there are no poles.

function bsmeasure(α::Function, n::Integer)
  (p, q, norms) = opuc(α, n, 1. + 0.*im)
  ϕ = x -> (1/sqrt(norms[n]))*polyval(p[n],x)
  return (θ -> 1./(τ*abs(ϕ(exp(im*θ)))^2))
end

end

