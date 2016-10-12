module Floquet

function floqdisc(f, λ, P)
  function rhs(t, y)
    return [y[2], f(t)*y[1] - λ*y[1]]
  end
  t1, y1 = ode45(rhs, [1.,0.], [0., P])
  t2, y2 = ode45(rhs, [0.,1.], [0., P])
  0.5*(y1[end][1] + y2[end][2])
end

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

function rmatRL(a::Number, x::Number)
  A = Array(Float64, 2, 2)
  A[1,1] = x - a
  A[1,2] = -1
  A[2,1] = 1
  A[2,2] = 0
  A
end

function tmatRL(a::Array, x::Number)
  A = eye(2)
  for i = 1:length(a)
    A *= rmatRL(a[i], x)
  end
  A
end

function discRL(pot::Array, x::Number)
  n = length(pot)
  A = SymTridiagonal(pot, ones(typeof(pot[1]),n-1)) - SymTridiagonal(fill(x, n), zeros(n-1))
  B = SymTridiagonal(pot[2:n-1], ones(typeof(pot[1]), n-3)) - SymTridiagonal(fill(x, n-2), zeros(n-3))
  return (det(A) - det(B))
end

#Incomplete, fix this
function rmatGC(α::Number, z::Number, j::Int64)
  A = Array(typeof(α), 2, 2)
  A[1,1] = z
  A[1,2] = -1 * conj(α) * z^j
  A[2,1] = -α * z^(-j)
  A[2,2] = z^-1
  A /=  #Need to fill in the correct norming constant here!
  A
end

function rmatUC(α::Number, z::Number)
  A = Array(typeof(α), 2, 2)
  A[1,1] = z
  A[1,2] = -1 * conj(α)
  A[2,1] = -α * z
  A[2,2] = 1
  A /= sqrt(1 - α*conj(α))
  A
end

function tmatUC(α::Array, z::Number)
  n = length(α)
  A = eye(2)
  for i = 1:n
    A *= rmatUC(α[i], z)
  end
  A
end

function tmatGC(α::Array, z::Number)
  A = tmatUC(α, z)
  A *= diagm([1, z^(-n)])
end

function discUC(α::Array, z::Number)
  n = length(α)
  if isodd(n)
    error("For technical reasons we assume even period")
  end
  return z^(-n/2) * trace(tmatUC(α, z))
end

end
