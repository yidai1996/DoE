# Batch
using OrdinaryDiffEq, LsqFit, XLSX, Plots, DifferentialEquations, NLsolve

function loadProcessData()
    global mu_maxE=1.7 #h^-1 from David's thesis(meeting slides from Prof.Lin) 1.7
    global mu_maxA=0.34 #h^-1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5478974/
    # global mu_maxA=0.091 #h^-1 https://www.sciencedirect.com/science/article/pii/S1369703X02001766
    # global mu_maxS=0.0504 #h^-1 https://onlinelibrary.wiley.com/doi/full/10.1002/cjce.22154
    global mu_maxS=0.5 #h^-1 https://onlinelibrary.wiley.com/doi/full/10.1002/cjce.22154
    global msE=0.1 # gsubstrate/gbiomass/h +-0.0008 h^-1 from David's thesis  substrate used for maintenence
    global msE_NH4=0.1 # gsubstrate/gbiomass/h
    global msA=0.31 # gsubstrate/gbiomass/h
    global msS=0.22 # gsubstrate/gbiomass/h
    global T0=303 #K
    global pH=7 # from Sofia medium
    global kdE=0.003 #h^-1 cell death rate from David's thesis
    global kdA=0.003 #h^-1
    global kdS=0.003 #h^-1
    global ksE=0.1 # gbiomass/L +-0.004 from David's thesis
    global ksA=0.1 # gbiomass/L https://www.sciencedirect.com/science/article/pii/S1369703X02001766
    global ksS=0.1 # gbiomass/L https://onlinelibrary.wiley.com/doi/full/10.1111/j.1529-8817.2005.04063.x
    global ksE_NH4=0.1
    # global ki=0.3 # my guessing
    global n=4.86
    global P_star=0.20
    global ysp_g=0.3 # gbiomass/gsubstrate from Minty 13 0.322
    global ysp_m=0.2 # gbiomass/gsubstrate from David's thesis, 0.409 at the first try
    global yspA=1.9 # https://microbialcellfactories.biomedcentral.com/track/pdf/10.1186/s12934-020-01362-9.pdf
    global yspS=2 #
    # global ysx=0.06 #http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=3 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    global ysxE=4 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    global ysxE_NH4=1.017
    global ysxA=0.160 # https://www.researchgate.net/figure/Growth-kinetics-of-Azotobacter-vinelandii-in-medium-before-and-after-optimization_tbl2_301753155
    global ysxS=0.17 # https://www.sciencedirect.com/science/article/pii/S0960852406004792
    # global D0= 0.68 # h^-1 initial dilusion rate
    global A0=0.01 # g/l initial substrate feeding concentration
    global E0=0.01 # g/L initial cell concentration 0.7g/L is from Figure 2.4 on Page 34 of David's thesis
    global S0=0.01 # g/L initial substrate concentration
    global N0=1 # g/L initial substrate concentration
    global C0=1 # g/L initial substrate concentration
    global P0=0.00 # g/L initial product concentration
    global Nin=0.2
    global Cin=0.2
    # global Xs=0.7 # g/L set point of cell concentration
    # global Ps=0.25 # g/L
    # global Ss=5.0
    global kO2=0.1 # guessing
    global M_O2=16 # g/mol
    global yox=1.017
    global mo=0.01 # guessing
    global kLaO=2.766*60 # 1/h # https://www.sciencedirect.com/science/article/pii/S0032959200002727
    global xO2_sat=7.5/16/1000 # mol/L https://www.waterboards.ca.gov/water_issues/programs/swamp/docs/cwt/guidance/3110en.pdf
    global O20=xO2_sat
    global tspan1=72 # h David's thesis

    println("Parameters Loaded!")

    xf=XLSX.readxlsx("G:\\My Drive\\Research\\DOE project\\Modeling\\Data fitting\\Monoculture with union media\\E.coli collected in October\\20221021_CN_Ec_OD.xlsx")
    sh1=xf["CN_OD"]
    global d1=convert(Array{Float64,2},sh1["A22:A28"])
    global E1=0.396*convert(Array{Float64,2},sh1["B22:B28"]) # C:N=5
    global E2=0.396*convert(Array{Float64,2},sh1["C22:C28"]) # C:N=15
    global E3=0.396*convert(Array{Float64,2},sh1["D22:D28"]) # C:N=36
    global E4=0.396*convert(Array{Float64,2},sh1["E22:E28"]) # C:N=60


    global X0=[E1[1] E2[1] E3[1] E4[1]] # g/L
    # global X0=[E11[2] E12[2] E21[1] E22[1]] # g/L
    global S0=[5 15 36 60] # g/L
    # global S0=[1 1 1 1] # g/L the assumed sucrose concentration at the second sample point
    global N0=[1 1 1 1] # g/L
    global PO=[.01 .01 .01 .01] # g/L
    println("Load Data Successfully!")
end

# global mumax=0.350
# global kd=0.006

function test1(tt,a)# Simple model: only sucrose limited
    # a=[mu_maxE ksE msE ysxE] kd=0。006
    # y=[E S P]

    # Without product term
    # a=[mu_maxE ksE msE ysxE] kd=0。006
    f(y,p,t)=[(a[1]*y[2]/(a[2]+y[2]) - kdE)*y[1], # X(E.coli)
                max(y[2],0)/y[2]*(-(a[1]*y[2]/(a[2]+y[2])/a[4] + a[3])*y[1])] # Sucrose
    # With product term
    # f(y,p,t)=[(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2]) - kdE)*y[1], # X(E.coli)
    #             max(y[2],0)/y[2]*(-(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])/a[4] + a[3])*y[1]), # Sucrose
    #             (max((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1],0)/((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2]) - kdE)*y[1])*a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])*ysp_g/a[4] + (1-max((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1],0)/(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1])*ysp_m*a[3])*y[1]] # Product
    # prob=ODEProblem(f,[X0[1],S0[1],PO[1]],(0,d1[end]))
    prob1=ODEProblem(f,[X0[1],S0[1]],(0,d1[end]))
    soln1=OrdinaryDiffEq.solve(prob1,Rosenbrock23())
    solns=zeros(size(tt)[1])
    for i=1:size(tt)[1]
        solns[i]=soln1(tt[i])[1]
    end
    println("solns=",solns)
    return solns
end

# function test2(tt,a) # Complicate model (two substrate limited: sucrose and ammonia)
#     # tt=times (assumed to be a vector), a = vector of parameters
#     # Complicate model (two substrate limited: sucrose and ammonia)
#     # a=[mu_maxE ksE ksE_NH4 msE msE_NH4] kd=0。006， mumax=0.230
#     # y=[E S N P]
#     f(y,p,t)=[(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3]) - kdE)*y[1], # X(E.coli)
#               max(y[2],0)/y[2]*(-(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])/ysxE + msE)*y[1]), # Sucrose
#               max(y[3],0)/y[3]*(- (mu_maxE*y[4]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])/ysxE_NH4 + msE_NH4)*y[1]), # Ammonia
#               (max((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1],0)/((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3]) - kdE)*y[1])*mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])*ysp_g/ysxE + (1-max((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1],0)/(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1])*ysp_m*msE)*y[1]] # Product
#
#
#     prob=ODEProblem(f,[X0[1],S0[1],PO[1]],(0,0,d1[end]))
#     solns=OrdinaryDiffEq.solve(prob,Rosenbrock23())
#     # TT=size(tt)[1]
#     # Tmax=17
#     # prob1=ODEProblem(f,[X0[1],S0[1]],(0.0,Tmax-T1[4]))
#     # soln1=OrdinaryDiffEq.solve(prob1,Rosenbrock23())
#     # prob2=ODEProblem(f,[X0[2],S0[2]],(0.0,Tmax-T1[4]))
#     # soln2=OrdinaryDiffEq.solve(prob2,Rosenbrock23())
#     # prob3=ODEProblem(f,[X0[3],S0[3]],(0.0,Tmax-T1[4]))
#     # soln3=OrdinaryDiffEq.solve(prob3,Rosenbrock23())
#     # prob4=ODEProblem(f,[X0[4],S0[4]],(0.0,Tmax-T1[4]))
#     # soln4=OrdinaryDiffEq.solve(prob4,Rosenbrock23())
#     # prob5=ODEProblem(f,[X0[5],S0[5]],(0.0,Tmax-T1[4]))
#     # soln5=OrdinaryDiffEq.solve(prob5,Rosenbrock23())
#     # solns=zeros(TT)
#     # for i=1:TT
#     #     if tt[i]<=17
#     #         solns[i]=soln5(tt[i]-T1[4])[1]/0.396 # give OD600 value
#     #     elseif tt[i]<=35
#     #         solns[i]=soln4(tt[i]-18-T1[4])[1]/0.396 # give OD600 value
#     #     elseif tt[i]<=53
#     #         solns[i]=soln3(tt[i]-36-T1[4])[1]/0.396 # give OD600 value
#     #     elseif tt[i]<=71
#     #         solns[i]=soln2(tt[i]-54-T1[4])[1]/0.396 # give OD600 value
#     #     else
#     #         solns[i]=soln1(tt[i]-72-T1[4])[1]/0.396# give OD600 value
#     #     end
#     #     # A=Array(soln)
#     #     # for j=1:size(solns)[1]
#     #     #     solns[j]=A[j]
#     #     # end
#     # end
#     return solns
# end


function ODEStep1(X,S,P,tspan,a) # Use one ODE solver to solve the whole system
    # f(y,p,t)=[(a[1]*y[2]/(a[2]+y[2]) - kdE)*y[1], # X(E.coli)
    #             max(y[2],0)/y[2]*(-(a[1]*y[2]/(a[2]+y[2])/a[4] + a[3])*y[1])] # Sucrose

    # No inhibition from isobutanol
    f(y,p,t)=[(a[1]*y[2]/(a[2]+y[2]) - kdE)*y[1], # X(E.coli)
                max(y[2],0)/y[2]*(-(a[1]*y[2]/(a[2]+y[2])/a[4] + a[3])*y[1]), # Sucrose
                (max((a[1]*y[2]/(a[2]+y[2])-kdE)*y[1],0)/((a[1]*y[2]/(a[2]+y[2]) - kdE)*y[1])*a[1]*y[2]/(a[2]+y[2])*ysp_g/a[4] + (1-max((a[1]*y[2]/(a[2]+y[2])-kdE)*y[1],0)/(a[1]*y[2]/(a[2]+y[2])-kdE)*y[1])*ysp_m*a[3])*y[1]] # Product
    # With inhibition
    # f(y,p,t)=[(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2]) - kdE)*y[1], # X(E.coli)
    #             max(y[2],0)/y[2]*(-(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])/a[4] + a[3])*y[1]), # Sucrose
    #             (max((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1],0)/((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2]) - kdE)*y[1])*a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])*ysp_g/a[4] + (1-max((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1],0)/(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1])*ysp_m*a[3])*y[1]] # Product
    prob1=ODEProblem(f,[X[1],S[1],P[1]],(0.0,tspan))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln1=DifferentialEquations.solve(prob1,Rosenbrock23())
    c1=soln1.t
    A1=Array(soln1)
    p1=plot(c1,A1[1,:],title="E.coli growth curve fitting",ylabel="E.coli Concentration (g/L)",xlabel="Time (hrs)",label="Fitting line of exper a",legend=:bottomright,framestyle=:box)
    scatter!(d1,E1,kind="scatter",label="Experimental data of a")
    display(p1)
    savefig("E.coli growth curve fitting.pdf")
    ps=plot(c1,A1[2,:],title="Sucrose consumed curve fitting",ylabel="Sucrose Concentration (g/L)",xlabel="Time (hrs)",label="Fitting line",legend=:bottomright,framestyle=:box)
    display(ps)
    savefig("Sucrose comsumption curve fitting.pdf")
    pp=plot(c1,A1[3,:],title="Isobutanol production without inhibition effection",ylabel="Isobutanol Concentration (g/L)",xlabel="Time (hrs)",label="Fitting line of exper a",legend=:bottomright,framestyle=:box)
    display(pp)
    savefig("Isobutanol production without inhibition effection curve fitting.pdf")
    return c1,A1[1,:],A1[2,:]
end

# function ODEStep2(X,S,N,P,tspan) # Use one ODE solver to solve the whole system
#     # Complicate
#     # f(y,p,t)=[(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3]) - kdE)*y[1], # X(E.coli)
#     #           max(y[2],0)/y[2]*(-(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])/ysxE + msE)*y[1]), # Sucrose
#     #           max(y[3],0)/y[3]*(- (mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE_NH4 + msE_NH4)*y[1]), # Ammonia
#     #           (max((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1],0)/((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3]) - kdE)*y[1])*mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])*ysp_g/ysxE + (1-max((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1],0)/(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1])*ysp_m*msE)*y[1]] # Product
#     # prob=ODEProblem(f,[X,S,N,P],(0.0,tspan))
#
#     # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
#     soln=DifferentialEquations.solve(prob,Rosenbrock23())
#     a=soln.t
#     A=Array(soln)
#     return a,A[1,:],A[2,:],A[3,:],A[4,:]
# end

function Datafit(Data,p,t)
    loadProcessData()
    # a=[mu_maxE ksE msE ysxE]
    # lower_p=[0.0,0.0,0.0,0.0,0.525]
    lower_p=[0.1,0.1,0.0,1]     #upper_p=[10,1000,1000,0.5,1]
    fit=curve_fit(test1,t,Data,p; lower=lower_p)#,upper=upper_p,maxIter=10000)
    Param=fit.param
    return Param
end

loadProcessData()
tt=vec(cat(d1,d1,dims=1))
CC=vec(cat(E4,E4,dims=1))
# lower_p=[0.0,0.0,0.0,0.0,0.525]
# p0=[.05,5,.2,.5,.56]
p0=[0.3,0.5,0.1,5]
Parameters=Datafit(CC,p0,tt)
ODEStep1(X0,S0,P0,d1[end],Parameters)
