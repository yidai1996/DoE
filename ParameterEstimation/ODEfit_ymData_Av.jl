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

    xf=XLSX.readxlsx("G:\\My Drive\\Research\\DOE project\\Modeling\\Data fitting\\Monoculture with union media\\Av\\Av_Mono_Av_Se_Flask.xlsx")
    sh1=xf["OD+pigments"]
    # sh2=xf["Iso_Prod_1"]
    # sh3=xf["OD_2"]
    # sh4=xf["Iso_Prod_2"]
    global d1=24*convert(Array{Float64,2},sh1["C18:C23"])
    global Av1=0.185*convert(Array{Float64,2},sh1["D18:D23"]) # g/L
    global Av2=0.185*convert(Array{Float64,2},sh1["E18:E23"])
    global N1=(14+4)*0.001*convert(Array{Float64,2},sh1["D27:D29"])
    global N2=(14+4)*0.001*convert(Array{Float64,2},sh1["E27:E29"])

    global X0=[Av1[1] Av2[1] N1[1] N2[1]] # g/L
    global S0=[5 5 5 5] # g/L
    # global S0=[1 1 1 1] # g/L the assumed sucrose concentration at the second sample point
    # global N0=[5 5 5 5] # g/L
    # global PO=[.01 .01 .01 .01] # g/L
    println("Load Data Successfully!")
end

# global mumax=0.350
# global kd=0.006

function test1(tt,a)# Simple model: only sucrose limited
    # a=[mu_maxA ksA msA ysxA kd] kd=0.006
    # y=[A S N]
    f(y,p,t)=[(a[1]*y[2]/(a[2]+y[2]) - a[5])*y[1],# X(Av)
    # f(y,p,t)=[(a[1]*y[2]/(a[2]+y[2]))*y[1],# X(Av)
         max(y[2],0)/y[2]*(- (a[1]*y[2]/(a[2]+y[2])/a[4]+a[3])*y[1])] # Sucrose
    # prob=ODEProblem(f,[X0[1],S0[1],PO[1]],(0,d1[end]))
    prob1=ODEProblem(f,[X0[1],S0[1]],(0,d1[end]))
    soln1=OrdinaryDiffEq.solve(prob1,Rosenbrock23())
    # prob2=ODEProblem(f,[X0[1],S0[1]],(0,d1[end]))
    # soln2=OrdinaryDiffEq.solve(prob2,Rosenbrock23())
    solns=zeros(size(tt)[1])
    for i=1:size(tt)[1]
        solns[i]=soln1(tt[i])[1]
    end
    println("solns=",solns)
    return solns
end

function ODEStep1(X,S,tspan,a) # Use one ODE solver to solve the whole system
    f(y,p,t)=[(a[1]*y[2]/(a[2]+y[2]) - a[5])*y[1], # X(Av)
                max(y[2],0)/y[2]*(-(a[1]*y[2]/(a[2]+y[2])/a[4] + a[3])*y[1])] # Sucrose
    prob1=ODEProblem(f,[X[1],S[1]],(0.0,tspan))
    prob2=ODEProblem(f,[X[2],S[2]],(0.0,tspan))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln1=DifferentialEquations.solve(prob1,Rosenbrock23())
    soln2=DifferentialEquations.solve(prob2,Rosenbrock23())
    c1=soln1.t
    A1=Array(soln1)
    c2=soln2.t
    A2=Array(soln2)
    p1=plot(c1,A1[1,:],title="A.v growth curve fitting",ylabel="A.v Concentration (g/L)",xlabel="Time (hrs)",label="Fitting line of exper a&b",legend=:bottomright,framestyle=:box)
    # p2=plot!(c2,A2[1,:],label="Fitting line of exper b")
    scatter!(d1,Av1,kind="scatter",label="Experimental data of a")
    scatter!(d1,Av2,kind="scatter",label="Experimental data of b")
    display(p1)
    savefig("A.v growth curve fitting.pdf")
    ps=plot(c1,A1[2,:],title="Sucrose consumed curve fitting",ylabel="Sucrose Concentration (g/L)",xlabel="Time (hrs)",label="Fitting line",legend=:bottomright,framestyle=:box)
    # display(ps)
    savefig("Sucrose comsumption curve fitting.pdf")
    return c1,c2,A1[1,:],A1[2,:],A2[1,:],A2[2,:]
end

function Datafit(Data,p,t)
    loadProcessData()
    CC=Data
    # a=[mu_maxE ksE msE ysxE]
    # lower_p=[0.0,0.0,0.0,0.0,0.525]
    lower_p=[0.0,0.0,0.0,0.0,0.0] 
    # lower_p=[0.0, 0.0, 0.0, 0.0]      #upper_p=[10,1000,1000,0.5,1]
    fit=curve_fit(test1,t,Data,p; lower=lower_p)#,upper=upper_p,maxIter=10000)
    Param=fit.param
    return Param
end

loadProcessData()
# tt=vec(cat(d1[2:5],d1[2:5],dims=1))
# CC=vec(cat(E11[2:5],E12[2:5],dims=1))
# tt=vec(cat(d1,d1,dims=1))
# CC=vec(cat(Av1,Av1,dims=1))
tt=vec(cat(d1,d1,dims=1))
CC=vec(cat(Av2,Av2,dims=1))
# lower_p=[0.0,0.0,0.0,0.0,0.525]
# p0=[.05,5,.2,.5,.56]
# p0=[0.3,0.5,0.1,5]
# p0=[0.108 0.787 0.52 2.27]
p0=[0.108, 0.52, 4, 2.27, 0.006] # should be a good starting point
# p0=[0.108, 0.52, 4, 2.27]
# ODEStep1(X0,S0,d1[end],p0)
p0=[ 0.10800000032712277, 0.5200000375563211, 4.000000000039041,2.2699642268716946,0.006000000288323959]
Parameters=Datafit(CC,p0,tt)
ODEStep1(X0,S0,d1[end],Parameters)
ODEStep1(X0,S0,d1[end],p0)
