using OrdinaryDiffEq, LsqFit, XLSX, Plots, DifferentialEquations, NLsolve

global mumax=0.230
global kd=0.006
# include(Monoculture_producer_Ecoli_ODEs_Batch_0223.jl)
xf=XLSX.readxlsx("G:\\My Drive\\Research\\DOE project\\Modeling\\Triculture\\roughly data fitting\\OD_Prod_Ec_2exp_datafit.xlsx")
sh1=xf["OD_1"]
sh2=xf["Iso_Prod_1"]
sh3=xf["OD_2"]
sh4=xf["Iso_Prod_2"]
d1=convert(Array{Float64,2},sh2["A16:A20"])
d2=convert(Array{Float64,2},sh4["B11:B14"])
E11=convert(Array{Float64,2},sh1["B12:B16"])
E12=convert(Array{Float64,2},sh1["C12:C16"])
E21=convert(Array{Float64,2},sh3["B11:B14"])
E22=convert(Array{Float64,2},sh3["C11:C14"])

P11=convert(Array{Float64,2},sh2["H16:H20"])
P12=convert(Array{Float64,2},sh2["I16:I20"])
P21=convert(Array{Float64,2},sh4["D11:D14"])
P22=convert(Array{Float64,2},sh4["G11:G14"])

global X0=[E11[1] E12[1] E21[1] E22[1]]*0.396 # g/L
global S0=[4 4 4 4] # g/L

# global mumax=0.350
# global kd=0.006
# function test(tt,a)
#     #tt=times (assumed to be a vector), a = vector of parameters
#     # a=[ks,ysx, ms] kd=0。006， mumax=0.230
#     f(y,p,t)=[(mumax*y[2]/(a[1]+y[2]) - kd)*y[1],
#          # -0.5*(tanh(100*y[2])+1)*(0.230*y[2]/(a[1]+y[2])/a[2] + a[3])*y[1]] # X,S
#          -max(y[2],0)/y[2]*(mumax*y[2]/(a[1]+y[2])/a[2] + a[3])*y[1]] # X,S
#     TT=size(tt)[1]
#     Tmax=17
# #    prob1=ODEProblem(f,[X0[1],S0[1]],(0.0,Tmax))
# #    soln1=OrdinaryDiffEq.solve(prob1,Rosenbrock23())
# #    prob2=ODEProblem(f,[X0[2],S0[2]],(0.0,Tmax))
# #    soln2=OrdinaryDiffEq.solve(prob2,Rosenbrock23())
#     prob3=ODEProblem(f,[X0[3],S0[3]],(0.0,Tmax))
#     soln3=OrdinaryDiffEq.solve(prob3,Rosenbrock23())
# #    prob4=ODEProblem(f,[X0[4],S0[4]],(0.0,Tmax))
# #    soln4=OrdinaryDiffEq.solve(prob4,Rosenbrock23())
#     prob5=ODEProblem(f,[X0[5],S0[5]],(0.0,Tmax))
#     soln5=OrdinaryDiffEq.solve(prob5,Rosenbrock23())
#     solns=zeros(TT)
#     for i=1:TT
#         if tt[i]<=17
#             solns[i]=soln3(tt[i])[1]/0.396 # give OD600 value
# #        elseif tt[i]<=35
# #            solns[i]=soln2(tt[i]-18)[1]/0.396 # give OD600 value
# #        elseif tt[i]<=53
# #            solns[i]=soln3(tt[i]-36)[1]/0.396 # give OD600 value
# #        elseif tt[i]<=71
# #            solns[i]=soln4(tt[i]-54)[1]/0.396 # give OD600 value
# #       else
# #            solns[i]=soln5(tt[i]-18)[1]/0.396# give OD600 value
#         end
#         # A=Array(soln)
#         # for j=1:size(solns)[1]
#         #     solns[j]=A[j]
#         # end
#     end
#     return solns
# end

function testL(tt,a) 
    #tt=times (assumed to be a vector), a = vector of parameters
    # a=[mu_max ks kO2] kd=0。006， mumax=0.230
    f(y,p,t)=[(a[4]*y[2]/(a[1]+y[2])*y[1])*(1-y[1]/a[5]) - kd*y[1],
         # -0.5*(tanh(100*y[2])+1)*(0.230*y[2]/(a[1]+y[2])/a[2] + a[3])*y[1]] # X,S
         -max(y[2],0)/y[2]*(mumax*y[2]/(a[1]+y[2])/a[2] + a[3])*y[1]] # X,S
    f(y,p,t)=[(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1]-D[1]*y[1], # X(E.coli)
              1/Vt*max(y[4],0)/y[4]*(-(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE + msE)*y[1] - ((mu_maxA*y[4]/(ksA+y[4]) - kdA)/ysxA+msA)*y[2] + yspS*(mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS*y[3]  - sum(D)*y[4]), # Sucrose
              1/Vt*max(y[5],0)/y[5]*(yspA*(mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]/ysxA - (mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE_NH4 + msE_NH4)*y[1] - ((mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS+msS)*y[3] - sum(D)*y[5]), # Ammonia
              1/Vt*(max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1])*mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])*ysp_g/ysxE + (1-max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1])*ysp_m*msE)*y[1] - sum(D)*y[6]] # Product

    TT=size(tt)[1]
    Tmax=17
    prob1=ODEProblem(f,[X0[1],S0[1]],(0.0,Tmax-T1[4]))
    soln1=OrdinaryDiffEq.solve(prob1,Rosenbrock23())
    prob2=ODEProblem(f,[X0[2],S0[2]],(0.0,Tmax-T1[4]))
    soln2=OrdinaryDiffEq.solve(prob2,Rosenbrock23())
    prob3=ODEProblem(f,[X0[3],S0[3]],(0.0,Tmax-T1[4]))
    soln3=OrdinaryDiffEq.solve(prob3,Rosenbrock23())
    prob4=ODEProblem(f,[X0[4],S0[4]],(0.0,Tmax-T1[4]))
    soln4=OrdinaryDiffEq.solve(prob4,Rosenbrock23())
    prob5=ODEProblem(f,[X0[5],S0[5]],(0.0,Tmax-T1[4]))
    soln5=OrdinaryDiffEq.solve(prob5,Rosenbrock23())
    solns=zeros(TT)
    for i=1:TT
        if tt[i]<=17
            solns[i]=soln5(tt[i]-T1[4])[1]/0.396 # give OD600 value
        elseif tt[i]<=35
            solns[i]=soln4(tt[i]-18-T1[4])[1]/0.396 # give OD600 value
        elseif tt[i]<=53
            solns[i]=soln3(tt[i]-36-T1[4])[1]/0.396 # give OD600 value
        elseif tt[i]<=71
            solns[i]=soln2(tt[i]-54-T1[4])[1]/0.396 # give OD600 value
        else
            solns[i]=soln1(tt[i]-72-T1[4])[1]/0.396# give OD600 value
        end
        # A=Array(soln)
        # for j=1:size(solns)[1]
        #     solns[j]=A[j]
        # end
    end
    return solns
end

function ODEStep(X,S,P,tspan) # Use one ODE solver to solve the whole system
    # f(y,p,t)=[(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2]) - kd)*y[1],
    f(y,p,t)=[(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1]-D[1]*y[1], # X(E.coli)
         (mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]-D[2]*y[2],# X(Av)
         (mu_maxS*y[5]/(ksS+y[5]) - kdS)*y[3]-D[3]*y[3],# X(Se)
         1/Vt*max(y[4],0)/y[4]*(-(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE + msE)*y[1] - ((mu_maxA*y[4]/(ksA+y[4]) - kdA)/ysxA+msA)*y[2] + yspS*(mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS*y[3]  - sum(D)*y[4]), # Sucrose
         1/Vt*max(y[5],0)/y[5]*(yspA*(mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]/ysxA - (mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE_NH4 + msE_NH4)*y[1] - ((mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS+msS)*y[3] - sum(D)*y[5]), # Ammonia
         1/Vt*(max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1])*mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])*ysp_g/ysxE + (1-max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1])*ysp_m*msE)*y[1] - sum(D)*y[6]] # Product
    prob=ODEProblem(f,[X,S,P],(0.0,tspan))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    return a,A[1,:],A[2,:],A[3,:]
end

function ODEStepL(X,S,P,tspan) # Use one ODE solver to solve the whole system
    # f(y,p,t)=[(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2]) - kd)*y[1],
    f(y,p,t)=[(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1]-D[1]*y[1], # X(E.coli)
         (mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]-D[2]*y[2],# X(Av)
         (mu_maxS*y[5]/(ksS+y[5]) - kdS)*y[3]-D[3]*y[3],# X(Se)
         1/Vt*max(y[4],0)/y[4]*(-(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE + msE)*y[1] - ((mu_maxA*y[4]/(ksA+y[4]) - kdA)/ysxA+msA)*y[2] + yspS*(mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS*y[3]  - sum(D)*y[4]), # Sucrose
         1/Vt*max(y[5],0)/y[5]*(yspA*(mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]/ysxA - (mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE_NH4 + msE_NH4)*y[1] - ((mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS+msS)*y[3] - sum(D)*y[5]), # Ammonia
         1/Vt*(max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1])*mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])*ysp_g/ysxE + (1-max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1])*ysp_m*msE)*y[1] - sum(D)*y[6]] # Product
    prob=ODEProblem(f,[X,S,P],(0.0,tspan))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    return a,A[1,:],A[2,:],A[3,:]
end

tt=vec(cat(T1[4:60],T1[4:60].+18,T1[4:60].+36,T1[4:60].+54,T1[4:60].+72,dims=1))
CC=vec(cat(C5[4:60],C4[4:60],C3[4:60],C2[4:60],C1[4:60],dims=1))
p0=[.05,5,.2,.5,.56] #a=[ks, ysx, ms,mumax,cap], best fit found by hand
lower_p=[0.0,0.0,0.0,0.0,0.525]
#upper_p=[10,1000,1000,0.5,1]
fit=curve_fit(testL,tt,CC,p0; lower=lower_p)#,upper=upper_p,maxIter=10000)
ks,ysx,ms,mumax,cap=fit.param

global tt1,X1,S1,P1=ODEStepL(C1[4]*0.396,2.0,0,T1[end]-T1[4])
global tt2,X2,S2,P2=ODEStepL(C2[4]*0.396,1.0,0,T1[end]-T1[4])
global tt3,X3,S3,P3=ODEStepL(C3[4]*0.396,0.75,0,T1[end]-T1[4])
global tt4,X4,S4,P4=ODEStepL(C4[4]*0.396,0.5,0,T1[end]-T1[4])
global tt5,X5,S5,P5=ODEStepL(C5[4]*0.396,0.25,0,T1[end]-T1[4])
plot(T1,[C1,C2,C3,C4,C5],xaxis="Time/h", yaxis="OD600",seriestype=:scatter,title="Experimental and model analysis of Ecoli.W (Sofia)",label="Experiment",legend=:topleft)
plot!([tt1.+T1[4],tt2.+T1[4],tt3.+T1[4],tt4.+T1[4],tt5.+T1[4]],[X1,X2,X3,X4,X5]/0.396,label="Model")

# tt=vec(T1)
# CC=vec(C2)
# p0=[1.1, 1.1, 1.1] #a=[ks, ysx, ms]
# lower_p=[0.0,0.0,0.0]
# fit=curve_fit(test,tt,CC,p0; lower=lower_p)
# ks,ysx,ms=fit.param
# global tt2,X2,S2,P2=ODEStep(C2[1]*0.396,1.0,0,tt[end])
# plot(T1,C2,xaxis="Time/h", yaxis="OD600",seriestype=:scatter,title="Experimental and model analysis of Ecoli.W (Sofia)",label="Experiment")
# plot!(tt2,X2/0.396,label="Model")
# plot(T3,C3,xaxis="Time/h", yaxis="OD600",seriestype=:scatter,title="Experimental and model analysis of Ecoli.W (Sofia)",label="Experiment")
# plot!(tt1,X1/0.396,label="Model")
# plot(T4,C4,xaxis="Time/h", yaxis="OD600",seriestype=:scatter,title="Experimental and model analysis of Ecoli.W (Sofia)",label="Experiment")
# plot!(tt1,X1/0.396,label="Model")
# plot(T5,C5,xaxis="Time/h", yaxis="OD600",seriestype=:scatter,title="Experimental and model analysis of Ecoli.W (Sofia)",label="Experiment")
# plot!(tt1,X1/0.396,label="Model")
