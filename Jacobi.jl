module Jacobi

importall ApproxFun
#importall SpecialMatrices
import Plots
import Roots
import Polynomials: Poly, polyval, roots

include("OPUC.jl")

function glmmat(F::Function, n::Int64, rank::Int64)
    entries = [F(n + i) for i=0:rank-1]
    M = Array{Float64}(rank, rank)
    for i = 1:rank
        for j = 1:rank
            M[i,j] = entries[1 + i-1 + j-1]
        end
    end
    M
end

function glmgmat(F::Function, n::Int64, rank::Int64)
    
end

function acoef(opF::Array, n::Integer)

end

function bcoef(opF::Array, n::Integer)

end

function invscatter(R::Function, poles::Array, normconst::Array)

end


