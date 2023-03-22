using LinearAlgebra
function JMatrix(mumax_Av, kd_Av, Ks_Av, Ys_Av, Y_NC_Av, C_N, X_Av, mumax_Se, kd_Se, Ks_Se, Ys_Se, Y_CN_Se, C_C, X_Se)
    J = zeros(4,4)
    mu_Av = mumax_Av*C_C/(Ks_Av+C_C)
    mu_Se = mumax_Se*C_N/(Ks_Se+C_N)
    J[1,1] = mu_Av - kd_Av
    J[1,2] = 0
    J[1,3] = mu_Av * Ks_Av /(Ks_Av+C_C)^2 *X_Av
    J[1,4] = 0
    J[2,1] = 0
    J[2,2] = mu_Se - kd_Se
    J[2,3] = 0
    J[2,4] = mu_Se * Ks_Se /(Ks_Se+C_N)^2 *X_Se
    J[3,1] = - mu_Av/Ys_Av
    J[3,2] = Y_CN_Se * mu_Se/Ys_Se
    J[3,3] = - X_Av / Ys_Av * mu_Av * Ks_Av /(Ks_Av+C_C)^2
    J[3,4] = Y_CN_Se * X_Se /Ys_Se *mu_Se * Ks_Se /(Ks_Se+C_N)^2
    J[4,1] = Y_NC_Av * mu_Av/Ys_Av
    J[4,2] = - mu_Se/Ys_Se
    J[4,3] = Y_NC_Av * X_Av /Ys_Av *mu_Av * Ks_Av /(Ks_Av+C_C)^2
    J[4,4] = - X_Se / Ys_Se * mu_Se * Ks_Se /(Ks_Se+C_N)^2

    return J
end

function StabilityDiagramDataGenerator(test)
    A = JMatrix(test[1], test[2], test[3], test[4], test[5], test[6], test[7], test[8], test[9], test[10], test[11], test[12], test[13], test[14])
    lambda = eigvals(A)
    println(lambda)
    vec = eigvecs(A)
    return lambda, vec
end

mumax_Av = 0.1
kd_Av = 0.1
Ks_Av = 0.3
Ys_Av = 1
Y_NC_Av = 1/3.5
# C_N = 2
# X_Av = 0
mumax_Se = 0.2
kd_Se = 0.1
Ks_Se = 0.01
Ys_Se = 5
Y_CN_Se = 3.5
# C_C = 0
# X_Se = 0
C_C = 2
C_N = 2
X_Av = 1
X_Se = X_Av*kd_Av/kd_Se*Ys_Se/Ys_Av*Y_NC_Av

test = [mumax_Av, kd_Av, Ks_Av, Ys_Av, Y_NC_Av, C_N, X_Av, mumax_Se, kd_Se, Ks_Se, Ys_Se, Y_CN_Se, C_C, X_Se]
lambda = StabilityDiagramDataGenerator(test)
include("Monoculture_producer_TriCulture_Chemostat_Monod.jl")
ParameterS = [mumax_Se, Ks_Se, Ys_Se, Y_CN_Se, 0]
ParameterA = [mumax_Av, Ks_Av, Ys_Av, Y_NC_Av, 0]
CocultureGrowth(ParameterS, ParameterA)

# TODO create the diagram for only change N0 and C0 (or X_Av and X_Se)