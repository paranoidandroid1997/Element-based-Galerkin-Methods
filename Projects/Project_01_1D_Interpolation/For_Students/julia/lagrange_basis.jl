function lagrange_basis(Np,Ns,ξ,xs)
    L = zeros(Np, Ns)
    dL = zeros(Np, Ns)
    for l = 1:Ns
        for i = 1:Np
            L[i, l] = 1.0
            for j = 1:Np
                prod = 1
                if (j != i)
                        L[i, l] *= (xs[l] - ξ[j])/(ξ[i] - ξ[j])
                    for k = 1:Np
                        if (k != i && k != j)
                            prod *= (xs[l] - ξ[k])/(ξ[i] - ξ[k])
                        end
                    end
                    dL[i,l] += prod/(ξ[i] - ξ[j])
                end
            end
        end
    end
    return L, dL
end #function