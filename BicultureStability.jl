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

function StabilityDiagramDataGenerator()
    A = JMatrix()
    lambda = eigvals(A)
    vec = eigvecs(A)
    return lambda, vec
end