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

function discriminant(pot::Array, z::Number)
  n = length(pot)
  truncpot = pot[2:n-1]
  A = SymTridiagonal(pot, ones(typeof(pot[1]),n-1)) - SymTridiagonal(fill(z, n), zeros(n-1))
  B = SymTridiagonal(truncpot, ones(typeof(pot[1]), n-3)) - SymTridiagonal(fill(z, n-2), zeros(n-3))
  return (det(A) - det(B))
end

function scatter()

end



end
