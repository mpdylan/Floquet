module Floquet

function permat(v::Array, h)
  n = 2*length(v)
  C = zeros(typeof(h), n, n)
  for i = 1:n-1 
    C[i, i+1] = 1 
  end
  C += transpose(C)
  for i = 1:div(n,2)
    C[i, i] = v[i]
    C[i + div(n,2), i + div(n,2)] = v[i]
  end
  C[1, n] = h
  C[n, 1] = 1/h
  return C
end

function floqdisc(v::Array, z::Number)
  n = length(v)
  truncv = v[2:n-1]
  A = SymTridiagonal(v, ones(typeof(v[1]),n-1)) - SymTridiagonal(fill(z, n), zeros(n-1))
  B = SymTridiagonal(truncv, ones(typeof(v[1]), n-3)) - SymTridiagonal(fill(z, n-2), zeros(n-2))
  return (det(A) - det(B))
end



end
