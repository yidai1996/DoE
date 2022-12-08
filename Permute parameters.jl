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
    # NA=[(mu[2+3]-mu[1+3])/mu[3+3], (ks[2+3]-ks[1+3])/ks[3+3], (ysx[2+3]-ysx[1+3])/ysx[3+3], (ysp[2+3]-ysp[1+3])/ysp[3+3], (ms[2+3]-ms[1+3])/ms[3+3]]
        
    original_listS=[mu[begin],ks[begin],ysx[begin],ysp[begin],ms[begin]]
    println(original_listS)
    original_listA=[mu[begin+3],ks[begin+3],ysx[begin+3],ysp[begin+3],ms[begin+3]]
    # parameter_list=zeros(Float64,(N[1],N[2],N[3],N[4],N[5],size(original_list)[1]))
    # for i1 in 0:N[1]
    #     for i2 in 0:N[2]
    #         for i3 in 0:N[3]
    #             for i4 in 0:N[4]
    #                 for i5 in 0:N[5]
    #                     parameter_list[i1,i2,i3,i4,i5,:]=original_list.+[i1*mu[end],i2*ks[end],i3*ysx[end],i4*ysp[end],i5*ms[end]]
    #                 end
    #             end
    #         end
    #     end
    # end

    parameter_listS = zeros(round(Int,prod(NS)),5)
    # println(parameter_listS)
    for i1 in 0:NS[1]-1
        for i2 in 0:NS[2]-1
            for i3 in 0:NS[3]-1
                for i4 in 0:NS[4]-1
                    for i5 in 1:NS[5]
                        parameter_listS[i1*NS[2]*NS[3]*NS[4]*NS[5]+i2*NS[3]*NS[4]*NS[5]+i3*NS[4]*NS[5]+i4*NS[5]+i5,:]=transpose(original_listS+[i1*mu[3],i2*ks[3],i3*ysx[3],i4*ysp[3],(i5-1)*ms[3]])
                        # parameter_listS[i1*NS[2]*NS[3]*NS[4]*NS[5]+i2*NS[3]*NS[4]*NS[5]+i3*NS[4]*NS[5]+i4*NS[5]+i5,:]=[0.02 0.02 0.02 0.02 0.02]
                    end
                end
            end
        end
    end
    # println(parameter_listS)
    # println(size(parameter_listS))
    num_permutationsS = prod(NS)
    println("number of permutation=",num_permutationsS)

    # Progress Bar
    for i in ProgressBar(1:num_permutationsS)
        CocultureGrowth(parameter_listS[i,:])
        
        # CocultureGrowth() will plot and store generated data.
    end

end
    