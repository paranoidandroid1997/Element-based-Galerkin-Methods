function element_matrices(ψ,dψ,Np,Nq,ωq,DFloat)
    dx = 
    M = zeros(Np, Np)

    for k = 1:Nq
        for i = 1:Np
            for j = 1:Np
                M[i,j] += ωq[k] * ψ[i,k] * ψ[j,k]
            end
        end
    end

    println(M)
    return M
end #function
