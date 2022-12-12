# Permute parameters of bicoculture
using Plots, DifferentialEquations, NLsolve, XLSX, DataFrames, ProgressBars, Profile, MathOptInterface, Printf
include("Monoculture_producer_TriCulture_Chemostat_Monod.jl")

function permute_parameters(mu,ks,ysx,ysp,ms)
    # mu=[lower_bound, upper_bound, # of interval]
    # Set up the permution parameters array
    # mu_list=collect(mu[begin]:mu[end]:mu[2])
    # ks_list=collect(ks[begin]:ks[end]:ks[2])
    # ysx_list=collect(ysx[begin]:ysx[end]:ysx[2])
    # ysp_list=collect(ysp[begin]:ysp[end]:ysp[2])
    # ms_list=collect(ms[begin]:ms[end]:ms[2])
    # N=[size(mu_list)[1], size(ks_list)[1], size(ysx_list)[1], size(ysp_list)[1], size(ms_list)[1]]
    NS=[round(Int,(mu[2]-mu[1])/mu[3])+1, round(Int,(ks[2]-ks[1])/ks[3])+1, round(Int,(ysx[2]-ysx[1])/ysx[3])+1, round(Int,(ysp[2]-ysp[1])/ysp[3])+1, round(Int,(ms[2]-ms[1])/ms[3])+1]
    println("NS=",NS)
    NA=[round(Int,(mu[2+3]-mu[1+3])/mu[3+3])+1, round(Int,(ks[2+3]-ks[1+3])/ks[3+3])+1, round(Int,(ysx[2+3]-ysx[1+3])/ysx[3+3])+1, round(Int,(ysp[2+3]-ysp[1+3])/ysp[3+3])+1, round(Int,(ms[2+3]-ms[1+3])/ms[3+3])+1]
        
    original_listS=[mu[begin],ks[begin],ysx[begin],ysp[begin],ms[begin]]
    println(original_listS)
    original_listA=[mu[begin+3],ks[begin+3],ysx[begin+3],ysp[begin+3],ms[begin+3]]
    println(original_listA)
    # parameter_list=zeros(Float64,(N[1],N[2],N[3],N[4],N[5],size(original_list)[1]))
    # for s1 in 0:N[1]
    #     for s2 in 0:N[2]
    #         for s3 in 0:N[3]
    #             for s4 in 0:N[4]
    #                 for s5 in 0:N[5]
    #                     parameter_list[s1,s2,s3,s4,s5,:]=original_list.+[s1*mu[end],s2*ks[end],s3*ysx[end],s4*ysp[end],s5*ms[end]]
    #                 end
    #             end
    #         end
    #     end
    # end

    parameter_listS = zeros(round(Int,prod(NS)),5)
    parameter_listA = zeros(round(Int,prod(NA)),5)
    # println(parameter_listS)
    for s1 in 0:NS[1]-1
        for s2 in 0:NS[2]-1
            for s3 in 0:NS[3]-1
                for s4 in 0:NS[4]-1
                    for s5 in 1:NS[5]
                        for a1 in 0:NA[1]-1
                            for a2 in 0:NA[2]-1
                                for a3 in 0:NA[3]-1
                                    for a4 in 0:NA[4]-1
                                        for a5 in 1:NA[5]
                                            parameter_listS[s1*NS[2]*NS[3]*NS[4]*NS[5]+s2*NS[3]*NS[4]*NS[5]+s3*NS[4]*NS[5]+s4*NS[5]+s5,:]=transpose(original_listS+[s1*mu[3],s2*ks[3],s3*ysx[3],s4*ysp[3],(s5-1)*ms[3]])
                                            parameter_listA[a1*NA[2]*NA[3]*NA[4]*NA[5]+a2*NA[3]*NA[4]*NA[5]+a3*NA[4]*NA[5]+a4*NA[5]+a5,:]=transpose(original_listA+[a1*mu[3+3],a2*ks[3+3],a3*ysx[3+3],a4*ysp[3+3],(a5-1)*ms[3+3]])
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    println(parameter_listA)
    # println(size(parameter_listS))
    num_permutationS = prod(NS)
    num_permutationA = prod(NA)
    # println("number of permutation=",num_permutationS)

    # Progress Bar
    for i in ProgressBar(1:num_permutationS)
        for j in ProgressBar(1:num_permutationA)
            CocultureGrowth(parameter_listS[i,:],parameter_listA[j,:])
        end        
        # CocultureGrowth() will plot and store generated data.
    end

end
    