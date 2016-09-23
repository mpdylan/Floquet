module Floquet

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
  truncpot = pot[2:n-1]
  A = SymTridiagonal(pot, ones(typeof(pot[1]),n-1)) - SymTridiagonal(fill(x, n), zeros(n-1))
  B = SymTridiagonal(truncpot, ones(typeof(pot[1]), n-3)) - SymTridiagonal(fill(x, n-2), zeros(n-3))
  return (det(A) - det(B))
end

function rmatUC(α::Complex, z::Complex)

  A = Array(typeof(α), 2, 2)
  A[1,1] = z
  A[1,2] = -1 * conj(α)
  A[2,1] = -α * z
  A[2,2] = 1
  A

end

function tmatUC(α::Array, z::Complex)

  n = length(α)
  A = eye(2)
  for i = 1:n
    A *= rmatUC(α[i], z)
  end
  A
end

function discUC(α::Array, z::Number)

  n = length(α)
  if isodd(n)
    error("For technical reasons we assume even period")
  end
  return z^(-p/2) * trace(tmatUC(α, z))
end

end
