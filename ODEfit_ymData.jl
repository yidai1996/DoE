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

    xf=XLSX.readxlsx("G:\\My Drive\\Research\\DOE project\\Modeling\\Triculture\\roughly data fitting (monoculture without union media)\\OD_Prod_Ec_2exp_datafit.xlsx")
    sh1=xf["OD_1"]
    sh2=xf["Iso_Prod_1"]
    sh3=xf["OD_2"]
    sh4=xf["Iso_Prod_2"]
    global d1=24*convert(Array{Float64,2},sh2["A16:A20"])
    global d2=24*convert(Array{Float64,2},sh4["B11:B14"])
    global E11=convert(Array{Float64,2},sh1["B12:B16"])
    global E12=convert(Array{Float64,2},sh1["C12:C16"])
    global E21=convert(Array{Float64,2},sh3["B11:B14"])
    global E22=convert(Array{Float64,2},sh3["C11:C14"])

    global P11=convert(Array{Float64,2},sh2["H16:H20"])
    global P12=convert(Array{Float64,2},sh2["I16:I20"])
    global P21=convert(Array{Float64,2},sh4["D11:D14"])
    global P22=convert(Array{Float64,2},sh4["G11:G14"])

    global X0=[E11[1] E12[1] E21[1] E22[1]]*0.396 # g/L
    global S0=[36 36 36 36] # g/L
    global N0=[5 5 5 5] # g/L
    global PO=[.01 .01 .01 .01] # g/L
    println("Load Data Successfully!")
end

# global mumax=0.350
# global kd=0.006

function test1(tt,a)# Simple model: only sucrose limited
    # a=[mu_maxE ksE msE ysxE] kd=0。006， mumax=0.230
    # y=[E S P]
    f(y,p,t)=[(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2]) - kdE)*y[1], # X(E.coli)
                max(y[2],0)/y[2]*(-(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])/a[4] + a[3])*y[1]), # Sucrose
                (max((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1],0)/((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2]) - kdE)*y[1])*a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])*ysp_g/a[4] + (1-max((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1],0)/(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1])*ysp_m*a[3])*y[1]] # Product
    prob=ODEProblem(f,[X0[1],S0[1],PO[1]],(0,d1[end]))
    soln=OrdinaryDiffEq.solve(prob,Rosenbrock23())
    # Tmax=17
    # prob1=ODEProblem(f,[X0[1],S0[1]],(0.0,Tmax-T1[4]))
    # soln1=OrdinaryDiffEq.solve(prob1,Rosenbrock23())
    # prob2=ODEProblem(f,[X0[2],S0[2]],(0.0,Tmax-T1[4]))
    # soln2=OrdinaryDiffEq.solve(prob2,Rosenbrock23())
    # prob3=ODEProblem(f,[X0[3],S0[3]],(0.0,Tmax-T1[4]))
    # soln3=OrdinaryDiffEq.solve(prob3,Rosenbrock23())
    # prob4=ODEProblem(f,[X0[4],S0[4]],(0.0,Tmax-T1[4]))
    # soln4=OrdinaryDiffEq.solve(prob4,Rosenbrock23())
    # prob5=ODEProblem(f,[X0[5],S0[5]],(0.0,Tmax-T1[4]))
    # soln5=OrdinaryDiffEq.solve(prob5,Rosenbrock23())
    solns=zeros(size(tt)[1])
    # for i=1:TT
    #     if tt[i]<=17
    #         solns[i]=soln5(tt[i]-T1[4])[1]/0.396 # give OD600 value
    #     elseif tt[i]<=35
    #         solns[i]=soln4(tt[i]-18-T1[4])[1]/0.396 # give OD600 value
    #     elseif tt[i]<=53
    #         solns[i]=soln3(tt[i]-36-T1[4])[1]/0.396 # give OD600 value
    #     elseif tt[i]<=71
    #         solns[i]=soln2(tt[i]-54-T1[4])[1]/0.396 # give OD600 value
    #     else
    #         solns[i]=soln1(tt[i]-72-T1[4])[1]/0.396# give OD600 value
    #     end
        A=Array(soln)/0.396
        for j=1:size(solns)[1]
            solns[j]=A[j]
        end
        println("solns=",solns)
    # end
    return solns
end

function test2(tt,a) # Complicate model (two substrate limited: sucrose and ammonia)
    # tt=times (assumed to be a vector), a = vector of parameters
    # Complicate model (two substrate limited: sucrose and ammonia)
    # a=[mu_maxE ksE ksE_NH4 msE msE_NH4] kd=0。006， mumax=0.230
    # y=[E S N P]
    f(y,p,t)=[(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3]) - kdE)*y[1], # X(E.coli)
              max(y[2],0)/y[2]*(-(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])/ysxE + msE)*y[1]), # Sucrose
              max(y[3],0)/y[3]*(- (mu_maxE*y[4]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])/ysxE_NH4 + msE_NH4)*y[1]), # Ammonia
              (max((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1],0)/((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3]) - kdE)*y[1])*mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])*ysp_g/ysxE + (1-max((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1],0)/(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1])*ysp_m*msE)*y[1]] # Product


    prob=ODEProblem(f,[X0[1],S0[1],PO[1]],(0,0,d1[end]))
    solns=OrdinaryDiffEq.solve(prob,Rosenbrock23())
    # TT=size(tt)[1]
    # Tmax=17
    # prob1=ODEProblem(f,[X0[1],S0[1]],(0.0,Tmax-T1[4]))
    # soln1=OrdinaryDiffEq.solve(prob1,Rosenbrock23())
    # prob2=ODEProblem(f,[X0[2],S0[2]],(0.0,Tmax-T1[4]))
    # soln2=OrdinaryDiffEq.solve(prob2,Rosenbrock23())
    # prob3=ODEProblem(f,[X0[3],S0[3]],(0.0,Tmax-T1[4]))
    # soln3=OrdinaryDiffEq.solve(prob3,Rosenbrock23())
    # prob4=ODEProblem(f,[X0[4],S0[4]],(0.0,Tmax-T1[4]))
    # soln4=OrdinaryDiffEq.solve(prob4,Rosenbrock23())
    # prob5=ODEProblem(f,[X0[5],S0[5]],(0.0,Tmax-T1[4]))
    # soln5=OrdinaryDiffEq.solve(prob5,Rosenbrock23())
    # solns=zeros(TT)
    # for i=1:TT
    #     if tt[i]<=17
    #         solns[i]=soln5(tt[i]-T1[4])[1]/0.396 # give OD600 value
    #     elseif tt[i]<=35
    #         solns[i]=soln4(tt[i]-18-T1[4])[1]/0.396 # give OD600 value
    #     elseif tt[i]<=53
    #         solns[i]=soln3(tt[i]-36-T1[4])[1]/0.396 # give OD600 value
    #     elseif tt[i]<=71
    #         solns[i]=soln2(tt[i]-54-T1[4])[1]/0.396 # give OD600 value
    #     else
    #         solns[i]=soln1(tt[i]-72-T1[4])[1]/0.396# give OD600 value
    #     end
    #     # A=Array(soln)
    #     # for j=1:size(solns)[1]
    #     #     solns[j]=A[j]
    #     # end
    # end
    return solns
end


function ODEStep1(X,S,P,tspan,a) # Use one ODE solver to solve the whole system
    f(y,p,t)=[(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2]) - kdE)*y[1], # X(E.coli)
                max(y[2],0)/y[2]*(-(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])/a[4] + a[3])*y[1]), # Sucrose
                (max((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1],0)/((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2]) - kdE)*y[1])*a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])*ysp_g/a[4] + (1-max((a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1],0)/(a[1]*y[2]*(max(1-y[3]/P_star,0.0))^n/(a[2]+y[2])-kdE)*y[1])*ysp_m*a[3])*y[1]] # Product
    prob=ODEProblem(f,[X,S,P],(0.0,tspan))

    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    c=soln.t
    A=Array(soln)
    p=plot(c,A[1,:])
    scatter!(d1,E11,kind="scatter")
    display(p)
    return c,A[1,:],A[2,:],A[3,:],A[4,:]
end

function testode(a,b,c,d)
    f(y,p,t)=[d[1]*y[1],
            exp(d[2]*y[2]),
            -y[3]+d[3]]
    prob=ODEProblem(f,[a,b,c],(0.0,5))
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    c=soln.t
    A=Array(soln)
    return c,A[1,:],A[2,:],A[3,:]
end

function ODEStep2(X,S,N,P,tspan) # Use one ODE solver to solve the whole system
    # Complicate
    # f(y,p,t)=[(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3]) - kdE)*y[1], # X(E.coli)
    #           max(y[2],0)/y[2]*(-(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])/ysxE + msE)*y[1]), # Sucrose
    #           max(y[3],0)/y[3]*(- (mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE_NH4 + msE_NH4)*y[1]), # Ammonia
    #           (max((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1],0)/((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3]) - kdE)*y[1])*mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])*ysp_g/ysxE + (1-max((mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1],0)/(mu_maxE*y[2]*(max(1-y[4]/P_star,0.0))^n/(ksE+y[2])*y[3]/(ksE_NH4+y[3])-kdE)*y[1])*ysp_m*msE)*y[1]] # Product
    # prob=ODEProblem(f,[X,S,N,P],(0.0,tspan))

    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    return a,A[1,:],A[2,:],A[3,:],A[4,:]
end

function Datafit()
    loadProcessData()
    CC=E11
    # a=[mu_maxE ksE msE]
    p0=[0.3,1,.1] #a=[ks, ysx, ms,mumax,cap], best fit found by hand
    lower_p=[0.0,0.0,0.0]     #upper_p=[10,1000,1000,0.5,1]
    fit=curve_fit(test1,d1,CC,p0; lower=lower_p)#,upper=upper_p,maxIter=10000)
    Param=fit.param
    return Param
end

p0=[0.5,2,.1,10]
ODEStep1(X0[1],S0[1],P0[1],d1[5],p0)

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
