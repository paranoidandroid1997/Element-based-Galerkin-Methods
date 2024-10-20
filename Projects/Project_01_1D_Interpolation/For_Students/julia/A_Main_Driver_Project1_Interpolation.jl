#=
-------------------------------------------------------------------------------------------------------------
This file runs the 1D Interpolation using 
The interpolation points used are the following:
ipoints=1: Lobatto
ipoints=2: Legendre
ipoints=3: Chebyshev
ipoints=4: Equi-spaced

This is part of Project 1 described in Algorithm 3.1 in F.X. Giraldo's Introduction to Element-based Galerkin Methods using 
Tensor-Product Bases: Analysis, Algorithms, and Applications.

Written by F.X. Giraldo on July 6, 2021.
Department of Applied Mathematics
Naval Postgraduate School
Monterey, CA 93943
-------------------------------------------------------------------------------------------------------------
=#

using Plots, LinearAlgebra, FastGaussQuadrature, LaTeXStrings

include("QuadraturePoints.jl")

#Some Constants
DFloat = Float64
Quadrature_type = "fxg"
#Quadrature_type = "julia"
Nmin = 1
Nmax = 64
Ns = 100
c = π / 2
Npoints = 4
machine_zero = eps(DFloat)

#Allocate Arrays
Narray = zeros(Int64, Nmax)
l1_norm_interpolation = zeros(DFloat, Nmax, Npoints)
l2_norm_interpolation = zeros(DFloat, Nmax, Npoints)
l8_norm_interpolation = zeros(DFloat, Nmax, Npoints)

function main()

    @show(DFloat, Quadrature_type)

    #Loop through polynomial orders
    for ipoints = 1:Npoints
        inop = 0
        for N = Nmin:Nmax
            Q = N
            Np = N + 1
            inop += 1
            Narray[inop] = N

            #Select Points
            if Quadrature_type == "fxg"
                if ipoints == 1
                    (ξ, ω) = QuadraturePoints.lobatto_gauss(Np) #Lobatto
                elseif ipoints == 2
                    (ξ, ω) = QuadraturePoints.legendre_gauss(Np) #Legendre
                elseif ipoints == 3
                    (ξ, ω) = QuadraturePoints.chebyshev_gauss(Np) #Chebyshev
                elseif ipoints == 4
                    (ξ, ω) = QuadraturePoints.equispaced_points(Np) #Equi-spaced
                end
            elseif Quadrature_type == "julia"
                if ipoints == 1
                    (ξ, ω) = gausslobatto(Np) #Lobatto
                elseif ipoints == 2
                    (ξ, ω) = gausslegendre(Np) #Legendre
                elseif ipoints == 3
                    (ξ, ω) = gausschebyshev(Np) #Chebyshev
                elseif ipoints == 4
                    (ξ, ω) = QuadraturePoints.equispaced_points(Np) #Equi-spaced
                end
            end

            #--------------------------------------------------#
            #Interpolation
            #--------------------------------------------------#
            #Compute Sample Space
            xs = zeros(Ns)
            xs = range(-1, length=Ns, stop=1)
            println(length(ξ))

            #=-------------------------------------------------
            Students Include the Lagrange Basis Function
            (ψ,dψ) = QuadraturePoints.lagrange_basis(Np,Ns,ξ,xs)
            ----------------------------------------------------=#
            (ψ, dψ) = QuadraturePoints.lagrange_basis(Np, Ns, ξ, xs)

            #Compute Expansion Coefficients
            q_coeff = zeros(DFloat, Np)
            for i = 1:Np
                x = ξ[i]
                q_coeff[i] = cos(c * x)
            end #i

            #Compute Nth Order Interpolant
            qn = zeros(DFloat, Ns)
            for i = 1:Ns
                qsum = 0
                for j = 1:Np
                    qsum = qsum + ψ[j, i] * q_coeff[j]
                end #j
                qn[i] = qsum
            end #i


            #Compute Exact Solution
            qe = zeros(DFloat, Ns)
            for i = 1:Ns
                x = xs[i]
                qe[i] = cos(c * x)
            end #i

            if (N == 8 && ipoints == 2)
                # Plot with dashed line and hollow dots
                p = plot(xs, qn,
                    label="interp",
                    title=L"Plot of $f(x) = \cos(\frac{\pi}{2}x)$",
                    xlabel="x",
                    ylabel="f(x)",
                    linestyle=:dash,           # Set line style to dashed
                    marker=:circle,            # Set marker to hollow circles
                    markersize=5,              # Adjust marker size
                    markerstrokewidth=2,       # Make the marker outline thicker for hollow effect
                    markercolor=:white,        # Hollow effect with white marker fill
                    linecolor=:blue)           # Set line color to blue
                
                plot!(xs, qe, linewidth=4, label="exact")

                # Optionally, add more features to the plot (like another line or more customization)
                #plot!(x, cos, label="cos(x)", linestyle=:dash)  # Adding cosine with a dashed line

                # Save the plot
                savefig(p, "../plots/cos-interp-1.png")
            end


            #Compute L1; L2; & L8 Norm
            l1_norm_interpolation[inop, ipoints] = norm(qn - qe, 1) / norm(qe, 1)
            l2_norm_interpolation[inop, ipoints] = norm(qn - qe, 2) / norm(qe, 2)
            l8_norm_interpolation[inop, ipoints] = norm(qn - qe, Inf) / norm(qe, Inf)
        end #N

    end #ipoints

    p1 = plot(Narray, l1_norm_interpolation, xlabel="N", ylabel="Error Norm", legend=true, lw=3, yaxis=:log, label=["Lobatto" "Legendre" "Chebyshev" "Equispaced"], title="L1 Interpolation Error")
    p2 = plot(Narray, l2_norm_interpolation, xlabel="N", ylabel="Error Norm", legend=true, lw=3, yaxis=:log, label=["Lobatto" "Legendre" "Chebyshev" "Equispaced"], title="L2 Interpolation Error")
    p3 = plot(Narray, l8_norm_interpolation, xlabel="N", ylabel="Error Norm", legend=true, lw=3, yaxis=:log, label=["Lobatto" "Legendre" "Chebyshev" "Equispaced"], title="L∞ Interpolation Error")
    #figure_1 = plot(p1,p2,p3,layout = (1,3))
    figure_1 = plot(p2)
    savefig(figure_1, "../plots/interp1.png")
    display(figure_1)


    #Plot Interpolation
    println("Done") #output

end

#----------------------------------#
# Run the main function
#----------------------------------#
main()
