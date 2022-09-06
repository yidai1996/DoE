# Simulation for Continuously Triculturing System
# Simplest version: the limited substrates are sucrose and ammonia for Av/Se respectively, while E.coli consider sucrose and oxygen transfer limitation
using Plots, JuMP, Ipopt, DifferentialEquations, NLsolve

function loadProcessData()
    global mu_maxE=1.7 #h^-1 from David's thesis(meeting slides from Prof.Lin) 1.7
    global mu_maxA=0.34 #h^-1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5478974/
    # global mu_maxA=0.091 #h^-1 https://www.sciencedirect.com/science/article/pii/S1369703X02001766
    # global mu_maxS=0.0504 #h^-1 https://onlinelibrary.wiley.com/doi/full/10.1002/cjce.22154
    global mu_maxS=0.5 #h^-1 https://onlinelibrary.wiley.com/doi/full/10.1002/cjce.22154
    global msE=0.1 # gsubstrate/gbiomass/h +-0.0008 h^-1 from David's thesis  substrate used for maintenence
    global msE_NH4=0.1 # gsubstrate/gbiomass/h
    global msA=0.1 # gsubstrate/gbiomass/h
    global msS=0.1 # gsubstrate/gbiomass/h
    global T0=303 #K
    global pH=7 # from Sofia medium
    global kdE=0.003 #h^-1 cell death rate from David's thesis
    global kdA=0.003 #h^-1
    global kdS=0.003 #h^-1
    global ksE=0.1 # gbiomass/L +-0.004 from David's thesis
    global ksA=0.31 # gbiomass/L https://www.sciencedirect.com/science/article/pii/S1369703X02001766
    global ksS=0.22 # gbiomass/L https://onlinelibrary.wiley.com/doi/full/10.1111/j.1529-8817.2005.04063.x
    global ksE_NH4=0.1
    # global ki=0.3 # my guessing
    global n=4.86
    global P_star=0.20
    global ysp_g=0.3 # gbiomass/gsubstrate from Minty 13 0.322
    global ysp_m=0.2 # gbiomass/gsubstrate from David's thesis, 0.409 at the first try
    global yspA=0.9 # https://microbialcellfactories.biomedcentral.com/track/pdf/10.1186/s12934-020-01362-9.pdf
    global yspS=2 #
    # global ysx=0.06 #http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=3 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    global ysxE=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
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
    global tspan1=72
    global tspan2=100
    global saiin=0.2  # L/min
    global saioutE=0.1 # David's thesis
    global saioutS=0.1
    global saioutA=0.1
    global V=0.5 #L
    global Vt=3*V
    global DE=saioutE/V # 0.2
    global DA=saioutA/V # 0.2
    global DS=saioutS/V # 0.2
    global D0=[DE DA DS]
    println("Parameters Loaded!")
end

function AllGrowth() # Continuous flow
    loadProcessData()
    global tt1,At1,St1,Ct1,Nt1=Startup(D0,A0,S0,N0,C0,tspan1)
    global tt,Et,At,St,Pt,Ct,Nt=Tripartite(D0,E0,At1[end],St1[end],Nt1[end],Ct1[end],P0,tspan2)
    # E,A,S mean E.coli, Av, and Se
    global dEdt=zeros(size(tt)[1])
    global dAdt=zeros(size(tt)[1])
    global dSdt=zeros(size(tt)[1])
    global dCdt1=zeros(size(tt)[1])
    global dCdt2=zeros(size(tt)[1])
    global dCdt3=zeros(size(tt)[1])
    global dNdt1=zeros(size(tt)[1])
    global dNdt2=zeros(size(tt)[1])
    global dNdt3=zeros(size(tt)[1])
    global dPdt=zeros(size(tt)[1])
    for i=1:size(tt)[1]
        dEdt[i]=(mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i]) - kdE)*Et[i]-saiout*Et[i]
    end
    plot(tt,dEdt,title="Growth rate of E.coli",xaxis="Time(hr)",yaxis="g/L/h",label=false)
    savefig("E.coliGrowthRate_test.png")

    for i=1:size(tt)[1]
        dAdt[i]=(mu_maxA*Ct[i]/(ksA+Ct[i]) - kdA)*At[i] - saiout*At[i]
    end
    plot(tt,dAdt,title="Growth rate of Av",xaxis="Time(hr)",yaxis="g/L/h",label=false)
    savefig("AvGrowthRate_test.png")

    for i=1:size(tt)[1]
        dSdt[i]=(mu_maxS*Nt[i]/(ksS+Nt[i]) - kdS)*St[i] - saiout*St[i]
    end
    plot(tt,dSdt,title="Growth rate of Se",xaxis="Time(hr)",yaxis="g/L/h",label=false)
    savefig("SeGrowthRate_test.png")

    for i=1:size(tt)[1]
        dCdt1[i]= - (mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i])/ysxE + msE)*Et[i]
        dCdt2[i]= - ((mu_maxA*Ct[i]/(ksA+Ct[i]) - kdA)/ysxA+msA)*At[i]
        dCdt3[i]= + yspS*(mu_maxS*Nt[i]/(ksS+Nt[i]) - kdS)/ysxS
    end
    plot(tt,dCdt1,title="Sucrose consumed/produced rate",xaxis="Time(hr)",yaxis="g/L/h",label="Consumed by E.coli")
    plot!(tt,dCdt2,xaxis="Time(hr)",yaxis="g/L/h",label="Consumed by Av")
    plot!(tt,dCdt3,xaxis="Time(hr)",yaxis="g/L/h",label="Produced by Se")
    savefig("Sucrose consumedandproduced rate_test.png")

    for i=1:size(tt)[1]
        dNdt1[i]= yspA*(mu_maxA*Ct[i]/(ksA+Ct[i]) - kdA)*At[i]/ysxA
        dNdt2[i]= - (mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i])/ysxE + msE)*Et[i]
        dNdt3[i]= - ((mu_maxS*Nt[i]/(ksS+Nt[i]) - kdS)/ysxS+msS)*St[i]
    end
    plot(tt,dNdt1,title="Ammonia consumed/produced rate",xaxis="Time(hr)",yaxis="g/L/h",label="Produced by Av")
    plot!(tt,dNdt2,xaxis="Time(hr)",yaxis="g/L/h",label="Consumed by E.coli")
    plot!(tt,dNdt3,xaxis="Time(hr)",yaxis="g/L/h",label="Consumed by Se")
    savefig("Ammonia consumedandproduced rate_test.png")

    for i=1:size(tt)[1]
        zz=(mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i]) - kdE)*Et[i]
        dPdt[i]=(max(zz,0)/zz*mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i])*ysp_g/ysxE + (1-max(zz,0)/zz)*ysp_m*msE)*Et[i]
    end
    plot(tt,dPdt,title="Isobutanol produced rate",xaxis="Time(hr)",yaxis="g/L/h",label=false)
    savefig("Isobutanol produced rate_test.png")

    plot(tt,Et/0.396,title="E.coli concentration profile",xaxis="Time(hr)",yaxis="OD600",label=false)
    savefig("Ecoli_test.png")
    plot(tt,At,title="Av concentration profile",xaxis="Time(hr)",yaxis="Av(g/L)",label=false)
    savefig("Av_test.png")
    plot(tt,St,title="Se concentration profile",xaxis="Time(hr)",yaxis="Se(g/L)",label=false)
    savefig("Se_test.png")
    plot(tt,Ct,title="Sucrose concentration profile",xaxis="Time(hr)",yaxis="Sucrose(g/L)",label=false)
    savefig("Sucrose_test.png")
    plot(tt,Nt,title="Ammonia concentration profile",xaxis="Time(hr)",yaxis="Ammonia(g/L)",label=false)
    savefig("Ammonia_test.png")
    plot(tt,Pt,title="Product concentration profile",xaxis="Time(hr)",yaxis="Isobutanol(g/L)",label=false)
    savefig("Isobutanol_test.png")
end

function Startup(D,A,S,N,C,tspan1)
    f(y,p,t)=[(mu_maxA*y[3]/(ksA+y[3]) - kdA)*y[1]-D[2]*y[1],# X(Av)
         (mu_maxS*y[4]/(ksS+y[4]) - kdS)*y[3]-D[3]*y[2],# X(Se)
         1/Vt*max(y[3],0)/y[3]*( - ((mu_maxA*y[3]/(ksA+y[3]) - kdA)/ysxA+msA)*y[1] + yspS*(mu_maxS*y[4]/(ksS+y[4]) - kdS)/ysxS*y[2] - (D[2]+D[3])*y[3]), # Sucrose
         1/Vt*max(y[4],0)/y[4]*(yspA*(mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]/ysxA  - ((mu_maxS*y[4]/(ksS+y[4]) - kdS)/ysxS+msS)*y[2] - (D[2]+D[3])*y[4])] # Ammonia
    prob=ODEProblem(f,[E,A,S,N,C,P],(0.0,tspan))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    # plot(a,A[1,:],title="Microbial concentration profile",xaxis="Time(hr)",yaxis="X(g/L)",label=false)
    # plot(a,A[2,:],title="Substrate concentration profile",xaxis="Time(hr)",yaxis="S(g/L)",label=false)
    # plot(a,A[3,:],title="Product concentration profile",xaxis="Time(hr)",yaxis="P(g/L)",label=false)
    return a,A[1,:],A[2,:],A[3,:],A[4,:]
    # At,St,Ct,Nt
end

function Tripartite(D,E,A,S,N,C,P,tspan2) # Use one ODE solver to solve the whole system
    # D: dilusion rate
    # E/A/S: e.coli/Ntrigen/Carbon fixer  concentration
    # The following odes haven't been modified yet
    f(y,p,t)=[(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1]-saiout*y[1], # X(E.coli)
         (mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]-saiout*y[2],# X(Av)
         (mu_maxS*y[5]/(ksS+y[5]) - kdS)*y[3]-saiout*y[3],# X(Se)
         max(y[4],0)/y[4]*(-(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE + msE)*y[1] - ((mu_maxA*y[4]/(ksA+y[4]) - kdA)/ysxA+msA)*y[2] + yspS*(mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS*y[3] + saiin*Cin - saiout*y[4]), # Sucrose
         max(y[5],0)/y[5]*(yspA*(mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]/ysxA - (mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE + msE)*y[1] - ((mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS+msS)*y[3] + saiin*Nin - saiout*y[5]), # Ammonia
         (max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1])*mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])*ysp_g/ysxE + (1-max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1])*ysp_m*msE)*y[1] - saiout*y[6]] # Product

    prob=ODEProblem(f,[E,A,S,N,C,P],(0.0,tspan))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    # plot(a,A[1,:],title="Microbial concentration profile",xaxis="Time(hr)",yaxis="X(g/L)",label=false)
    # plot(a,A[2,:],title="Substrate concentration profile",xaxis="Time(hr)",yaxis="S(g/L)",label=false)
    # plot(a,A[3,:],title="Product concentration profile",xaxis="Time(hr)",yaxis="P(g/L)",label=false)
    return a,A[1,:],A[2,:],A[3,:],A[6,:],A[4,:],A[5,:]
    # Et,At,St,Pt,Ct,Nt
end
