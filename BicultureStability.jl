using LinearAlgebra
function JMatrix4D(mumax_Av, kd_Av, Ks_Av, Ys_Av, Y_NC_Av, C_N, X_Av, mumax_Se, kd_Se, Ks_Se, Ys_Se, Y_CN_Se, C_C, X_Se)
    J = zeros(4,4)
    # Only Monod model
    mu_Av = mumax_Av*C_C/(Ks_Av+C_C)
    mu_Se = mumax_Se*C_N/(Ks_Se+C_N)
    J[1,1] = mu_Av - kd_Av
    J[1,2] = 0
    J[1,3] = mu_Av * Ks_Av / (Ks_Av+C_C)^2 * X_Av
    J[1,4] = 0
    J[2,1] = 0
    J[2,2] = mu_Se - kd_Se
    J[2,3] = 0
    J[2,4] = mu_Se * Ks_Se /(Ks_Se+C_N)^2 * X_Se
    J[3,1] = - mu_Av / Ys_Av
    J[3,2] = Y_CN_Se * mu_Se / Ys_Se
    J[3,3] = - X_Av / Ys_Av * mumax_Av * Ks_Av /(Ks_Av+C_C)^2
    J[3,4] = Y_CN_Se * X_Se / Ys_Se *mumax_Se * Ks_Se /(Ks_Se+C_N)^2
    J[4,1] = Y_NC_Av * mu_Av / Ys_Av
    J[4,2] = - mu_Se / Ys_Se
    J[4,3] = Y_NC_Av * X_Av / Ys_Av *mumax_Av * Ks_Av /(Ks_Av+C_C)^2
    J[4,4] = - X_Se / Ys_Se * mumax_Se * Ks_Se /(Ks_Se+C_N)^2

    # Monod + L-P model
    # mu_Av = mumax_Av*C_C/(Ks_Av+C_C)
    # mu_Se = mumax_Se*C_N/(Ks_Se+C_N)
    # J[1,1] = mu_Av - kd_Av
    # J[1,2] = 0
    # J[1,3] = mu_Av * Ks_Av / (Ks_Av+C_C)^2 * X_Av
    # J[1,4] = 0
    # J[2,1] = 0
    # J[2,2] = mu_Se - kd_Se
    # J[2,3] = 0
    # J[2,4] = mu_Se * Ks_Se /(Ks_Se+C_N)^2 * X_Se
    # J[3,1] = - mu_Av / Ys_Av
    # J[3,2] = Y_CN_Se * (mu_Se - kd_Se) + beta_C
    # J[3,3] = - X_Av / Ys_Av * mumax_Av * Ks_Av / (Ks_Av + C_C)^2
    # J[3,4] = Y_CN_Se * X_Se *mumax_Se * Ks_Se / (Ks_Se + C_N)^2
    # J[4,1] = Y_NC_Av * (mu_Av -kd_Av) + beta_N
    # J[4,2] = - mu_Se / Ys_Se
    # J[4,3] = Y_NC_Av * X_Av *mumax_Av * Ks_Av / (Ks_Av + C_C)^2
    # J[4,4] = - X_Se / Ys_Se * mumax_Se * Ks_Se / (Ks_Se + C_N)^2
    return J
end

function StabilityDiagramDataGenerator4D(test)
    A = JMatrix4D(test[1], test[2], test[3], test[4], test[5], test[6], test[7], test[8], test[9], test[10], test[11], test[12], test[13], test[14])
    lambda = eigvals(A) # X_Av, X_Se, C_C, C_N
    # println(lambda)
    vec = eigvecs(A)
    # print(size(vec))
    return lambda, vec, A
end

# function JMatrix2D(mumax_Av, kd_Av, Ks_Av, Ys_Av, Y_NC_Av, C_N, X_Av, mumax_Se, kd_Se, Ks_Se, Ys_Se, Y_CN_Se, C_C, X_Se)
#     J = zeros(2,2)
#     mu_Av = mumax_Av*C_C/(Ks_Av+C_C)
#     mu_Se = mumax_Se*C_N/(Ks_Se+C_N)
#     J[1,1] = mu_Av - kd_Av
#     J[1,2] = 0
#     J[2,1] = 0
#     J[2,2] = mu_Se - kd_Se

#     # J[1,1] = - X_Av / Ys_Av * mumax_Av * Ks_Av /(Ks_Av+C_C)^2
#     # J[1,2] = Y_CN_Se * X_Se / Ys_Se *mumax_Se * Ks_Se /(Ks_Se+C_N)^2
#     # J[2,1] = Y_NC_Av * X_Av / Ys_Av *mumax_Av * Ks_Av /(Ks_Av+C_C)^2
#     # J[2,2] = - X_Se / Ys_Se * mumax_Se * Ks_Se /(Ks_Se+C_N)^2
#     # println(J)

#     return J
# end

# function StabilityDiagramDataGenerator2D(test)
#     A = JMatrix2D(test[1], test[2], test[3], test[4], test[5], test[6], test[7], test[8], test[9], test[10], test[11], test[12], test[13], test[14])
#     lambda = eigvals(A) # X_Av, X_Se, C_C, C_N
#     # println(lambda)
#     vec = eigvecs(A)
#     # print(size(vec))
#     return lambda, vec
# end





mumax_Av = 0.5
kd_Av = 0.05*mumax_Av # 5% of mumax_Av
# Ks_Av = 5*17.03*10^(-3)  # 5mM of sucrose
Ks_Av = 5*342.3*10^(-3) # 5mM of sucrose
Ys_Av = 1
Y_NC_Av = 2 # alpha_N
mumax_Se = 0.1
kd_Se = 0.05*mumax_Se
Ks_Se = 5*17.03*10^(-3) # 5mM of ammonia
# Ks_Se = 0.02 # 5mM of ammonia
# Ks_Se = 0.344 # 1.5 1.3 1.1 0.9 0.7 0.5 0.3 0.1
Ys_Se = 1  
Y_CN_Se = 0.5 # alpha_C
beta_C = 0 # 0 # 0.5
beta_N = 0 # 0 # 0.5

# ParameterS = [mumax_Se, Ks_Se, Ys_Se, Y_CN_Se, 0];
# ParameterA = [mumax_Av, Ks_Av, Ys_Av, Y_NC_Av, 0];
# C_N = 2
# C_C = 0
# CocultureGrowth(ParameterS, ParameterA, C_N, C_C; filename = "C_$(@sprintf("%.2f",C_C))_N__$(@sprintf("%.2f",C_N))")

# X_Ax = X_Se = 0
# C_N_list = collect(0:0.1:4);
C_C_list = LinRange(0,10,11);
C_N_list = 0.1*C_C_list;
X_Se = 0
X_Av = 0

Stability = zeros(length(C_C_list)*length(C_N_list),3);
println(length(C_C_list)*length(C_N_list))
include("Monoculture_producer_TriCulture_Chemostat_Monod.jl")
function StabilityAnalysis2D_Nutrient(C_N_list, C_C_list)
    eigenvalue = zeros(length(C_C_list)*length(C_N_list), 4);
    eigenvector = zeros(length(C_C_list)*length(C_N_list), 4, 4);
    # 2D case (biomass are not variables)
    # eigenvalue = zeros(length(C_C_list)*length(C_N_list), 2);
    # eigenvector = zeros(length(C_C_list)*length(C_N_list), 2, 2);
    count = 0
    println("start permutating C_N and C_C")
    for i in eachindex(C_N_list)
        C_N = C_N_list[i] 
        for j in eachindex(C_C_list)
            C_C = C_C_list[j]
            count += 1
            println("count=",count)
            test = [mumax_Av, kd_Av, Ks_Av, Ys_Av, Y_NC_Av, C_N, X_Av, mumax_Se, kd_Se, Ks_Se, Ys_Se, Y_CN_Se, C_C, X_Se]
            eigenvalue[count, :], eigenvector[count,:,:] = StabilityDiagramDataGenerator4D(test)
            # eigenvalue[count, :], eigenvector[count,:,:] = StabilityDiagramDataGenerator2D(test)
            # println("size of eigenvector = ",size(vec))
            # println(eigenvector[count,:,:])
            
            # Classify stability
            # 1 = stable 2 = Lyapunov stable 3 = unstable
            if sum(x->x>0, eigenvalue[count,:]) > 0
                Stability[count,:] = [C_C C_N 3]
            elseif sum(x->x==0, eigenvalue[count,:]) > 0
                Stability[count,:] = [C_C C_N 2]
            else Stability[count,:] = [C_C C_N 1]
            end

            # ParameterS = [mumax_Se, Ks_Se, Ys_Se, Y_CN_Se, 0.5]; # the last element is ms
            # ParameterA = [mumax_Av, Ks_Av, Ys_Av, Y_NC_Av, 0.5];
            ParameterS = [mumax_Se, Ks_Se, Ys_Se, Y_CN_Se, beta_C];
            ParameterA = [mumax_Av, Ks_Av, Ys_Av, Y_NC_Av, beta_N];
            CocultureGrowth(ParameterS, ParameterA, C_N, C_C; filename = "C_$(@sprintf("%.2f",C_C))_N__$(@sprintf("%.2f",C_N))")
        end
    end

    return eigenvalue, eigenvector, Stability
end

A,B,C = StabilityAnalysis2D_Nutrient(C_N_list, C_C_list)
using DataFrames, CSV
df = DataFrame(Sucrose = C[:,1], Ammonia = C[:,2], Stability = C[:,3], Eigenvalue1 = A[:,1], Eigenvalue2 = A[:,2], Eigenvalue3 = A[:,3], Eigenvalue4 = A[:,4], Eigenvector1_1 = B[:,1,1], Eigenvector1_2 = B[:,1,2], Eigenvector1_3 = B[:,1,3], Eigenvector1_4 = B[:,1,4], 
    Eigenvector2_1 = B[:,2,1], Eigenvector2_2 = B[:,2,2], Eigenvector2_3 = B[:,2,3], Eigenvector2_4 = B[:,2,4], Eigenvector3_1 = B[:,3,1], Eigenvector3_2 = B[:,3,2], Eigenvector3_3 = B[:,3,3], Eigenvector3_4 = B[:,3,4], Eigenvector4_1 = B[:,4,1], Eigenvector4_2 = B[:,4,2], Eigenvector4_3 = B[:,4,3], Eigenvector4_4 = B[:,4,4])
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\doe22_tspan_200\\StabilityAnalysisResults\\LinearStabilityRegionsEigenValueVectors.csv",df)

index_unstable = findall(x-> x==3, C[:,3])
index_lyapunov_stable = findall(x-> x==2, C[:,3])
index_stable = findall(x-> x==1, C[:,3])

show(stdout, "text/plain", A[index_stable,:])
show(stdout, "text/plain", A[index_lyapunov_stable,:])
show(stdout, "text/plain", A[index_unstable,:])

num_unstable = length(index_unstable)
num_lyapunov_stable = length(index_lyapunov_stable)
num_stable = length(index_stable)

Real_Unstable = zeros(num_unstable,2)
Real_lyapunov_stable = zeros(num_lyapunov_stable,2)
Real_stable = zeros(num_stable,2)

Real_Unstable = C[index_unstable,1:2]
Real_lyapunov_stable = C[index_lyapunov_stable,1:2]
Real_stable = C[index_stable,1:2]

show(stdout, "text/plain", C[index_stable,:])
show(stdout, "text/plain", C[index_lyapunov_stable,:])
show(stdout, "text/plain", C[index_unstable,:])

New_data_0 = C[index_unstable,1:2]
New_data_1 = C[index_lyapunov_stable,1:2]
New_data_2 = C[index_stable,1:2]

using Plots
scatter(New_data_0[:,1], New_data_0[:,2], xlabel = "C(Sucrose)", ylabel = "C(Ammonia)", label = "Unstable", markersize = 4, size = [1000,1000], markercolor = "Deep Sky Blue", markerstrokewidth = 0, legendfontsize=18, xtickfont = 18, ytickfont = 18, xguidefont = 18, yguidefont = 18)
scatter!(New_data_1[:,1], New_data_1[:,2], label = "Lyapunov Stable", markersize = 4,  markercolor = "Light Pink", markerstrokewidth = 0)
scatter!(New_data_2[:,1], New_data_2[:,2], label = "Stable", markersize = 4, markercolor = "Dark Magenta", markerstrokewidth = 0)
Plots.savefig("C:\\Users\\yid\\TemporaryResearchDataStorage\\doe21_tspan_200\\StabilityAnalysisResults\\LinearAnaylsis.pdf")

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
Plots.plot(rectangle(100, 100, 0, 0), xlabel="C(Sucrose)", ylabel="C(Ammonia)", label="Unstable")
# Plots.plot(rectangle(100, 100, 0, 0), xlabel="C(Sucrose)", ylabel="C(Ammonia)", label="Lyapunov")
Plots.plot!(Real_lyapunov_stable[:,1], Real_lyapunov_stable[:,2],linewidth=10, label="Lyapunov Stable")

using DataFrames
df = DataFrame([Float64[],Float64[],[]], ["C_C","C_N","Stability_Categories"])
for i in 1:length(C_C_list)*length(C_N_list)
    push!(df, C[i,:])
end
print(df)

# using Plots
# Plots.plot(C[index_stable,1], C[index_stable,2], xlabel="Sucrose", ylabel="Ammonia", label="Stable", seriestype=:scatter, markershape = :none, markerstrokewidth = 0, markercolor = :green)
# Plots.plot(C[index_lyapunov_stable,1], C[index_lyapunov_stable,2], label="Lyapunov Stable", seriestype=:scatter, markershape = :none, markerstrokewidth = 0, markercolor = :blue)
# Plots.plot(C[index_unstable,1], C[index_unstable,2], label="Stable", seriestype=:scatter, markershape = :none, markerstrokewidth = 0, markercolor = :red)

# using PlotlyJS
# PlotlyJS.plot(
#     df, x=:X_Av, y=:X_Se, color=:Stability_Categories,
#     mode="markers"
# )

# show(stdout,"text/plain", C)

# test if jacobian calculation is correct
using ForwardDiff
import ForwardDiff.jacobian
f(Av,Se,C,N) = [(mumax_Av*C/(Ks_Av+C) - kd_Av)*Av, 
                (mumax_Se*N/(Ks_Se+N) - kd_Se)*Se, 
                Y_CN_Se*(mumax_Se*N/(Ks_Se+N) - kd_Se)*Se + beta_C*Se - mumax_Av*C/(Ks_Av+C)*Av/Ys_Av,
                Y_NC_Av*(mumax_Av*C/(Ks_Av+C) - kd_Av)*Av + beta_N*Av - mumax_Se/Ys_Se*N/(Ks_Se+N)*Se]
x0 = [0.0, 0.0, 0.004715, 0.111833];
J0 = jacobian(x -> f(x...) ,x0)
# show(stdout, "text/plain", Real_lyapunov_stable)
