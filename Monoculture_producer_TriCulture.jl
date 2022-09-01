# Simulation for Continuously Triculturing System
# Simplest version: the limited substrates are sucrose and ammonia for Av/Se respectively, while E.coli consider sucrose and oxygen transfer limitation
using Plots, JuMP, Ipopt, DifferentialEquations, NLsolve

function loadProcessData()
    global mu_maxE=1.7 #h^-1 from David's thesis(meeting slides from Prof.Lin) 1.7
    global mu_maxA=1.7 #h^-1
    global mu_maxS=1.7 #h^-1
    global msE=1 # gsubstrate/gbiomass/h +-0.0008 h^-1 from David's thesis  substrate used for maintenence
    global msA=1 # gsubstrate/gbiomass/h
    global msS=1 # gsubstrate/gbiomass/h
    global T0=303 #K
    global pH=7 # from Sofia medium
    global kdE=0.003 #h^-1 cell death rate from David's thesis
    global kdA=0.003 #h^-1
    global kdS=0.003 #h^-1
    global ksE=0.1 # gbiomass/L +-0.004 from David's thesis
    global ksA=0.1 # gbiomass/L
    global ksS=0.1 # gbiomass/L
    # global ki=0.3 # my guessing
    global n=4.86
    global P_star=0.20
    global ysp_g=0.3 # gbiomass/gsubstrate from Minty 13 0.322
    global ysp_m=0.2 # gbiomass/gsubstrate from David's thesis, 0.409 at the first try
    global yspA=0.2 #
    global yspS=0.2 #
    # global ysx=0.06 #http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=3 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    global ysxE=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    global ysxA=1.017 #
    global ysxS=1.017 #
    # global D0= 0.68 # h^-1 initial dilusion rate
    global A0=0.01 # g/l initial substrate feeding concentration
    global E0=0.01 # g/L initial cell concentration 0.7g/L is from Figure 2.4 on Page 34 of David's thesis
    global S0=0.01 # g/L initial substrate concentration
    global N0=1 # g/L initial substrate concentration
    global C0=1 # g/L initial substrate concentration
    global P0=0.00 # g/L initial product concentration
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
    global tspan=100
    global saiin=0.1  # L/min
    global saiout=0.1
    global V=5 #L
    global D0=saiout/V
    println("Parameters Loaded!")
end

function AllGrowth() # Continuous flow
    loadProcessData()
    global tt,Et,At,St,Pt,Ct,Nt=Tripartite(D0,E0,A0,S0,N0,C0,P0,tspan)
    # E,A,S mean E.coli, Av, and Se
    # global dxdt=zeros(size(tt1)[1])
    # global dPdt=zeros(size(tt1)[1])
    # global dSdt=zeros(size(tt1)[1])
    # for i=1:size(tt1)[1]
    #     dxdt[i]=(mu_max*S1[i]*(max(1-P1[i]/P_star,0.0))^n/(ks+S1[i]) - kd)*X1[i]
    # end
    # for i=1:size(tt1)[1]
    #     dSdt[i]=- (mu_max*S1[i]*(max(1-P1[i]/P_star,0.0))^n/(ks+S1[i])/ysx + ms)*X1[i]
    # end
    # for i=1:size(tt1)[1]
    #     zz=(mu_max*S1[i]*(max(1-P1[i]/P_star,0.0))^n/(ks+S1[i])-kd)*X1[i]
    #     dPdt[i]=(max(zz,0)/zz*mu_max*S1[i]*(max(1-P1[i]/P_star,0.0))^n/(ks+S1[i])*ysp_g/ysx + (1-max(zz,0)/zz)*ysp_m*ms)*X1[i]
    # end
    plot(tt,Et/0.396,title="E.coli concentration profile",xaxis="Time(hr)",yaxis="OD600",label=false)
    savefig("Ecoli_test.png")
    # plot(tt1,dxdt,title="Microbial growth rate",xaxis="Time(hr)",yaxis="dXdt(g/L/h)",label=false)
    plot(tt,At,title="Av concentration profile",xaxis="Time(hr)",yaxis="Av(g/L)",label=false)
    savefig("Av_test.png")
    plot(tt,St,title="Se concentration profile",xaxis="Time(hr)",yaxis="Se(g/L)",label=false)
    savefig("Se_test.png")
    plot(tt,Ct,title="Sucrose concentration profile",xaxis="Time(hr)",yaxis="Sucrose(g/L)",label=false)
    savefig("Sucrose_test.png")
    plot(tt,Nt,title="Ammonia concentration profile",xaxis="Time(hr)",yaxis="Ammonia(g/L)",label=false)
    savefig("Ammonia_test.png")
    # plot(tt1,dSdt,title="Substrate change rate",xaxis="Time(hr)",yaxis="dSdt(g/L/h)",label=false)
    plot(tt,Pt,title="Product concentration profile",xaxis="Time(hr)",yaxis="Isobutanol(g/L)",label=false)
    savefig("Isobutanol_test.png")
    # plot(tt1,dPdt,title="Product rate",xaxis="Time(hr)",yaxis="dPdt(g/L/h)",label=false)
end

function Tripartite(D,E,A,S,N,C,P,tspan) # Use one ODE solver to solve the whole system
    # D: dilusion rate
    # E/A/S: e.coli/Ntrigen/Carbon fixer  concentration
    # The following odes haven't been modified yet
    f(y,p,t)=[(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*M_O2*O20/(kO2+M_O2*O20) - kdE)*y[1], # X(E.coli)
         (mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2],# X(Av)
         (mu_maxS*y[5]/(ksS+y[5]) - kdS)*y[3],# X(Se)
         max(y[4],0)/y[4]*(-(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])/ysxE + msE)*y[1] - ((mu_maxA*y[4]/(ksA+y[4]) - kdA)/ysxA+msA)*y[2]#=+saiin/V*Sin=#-saiout/V + yspS*(mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS), # Sucrose
         max(y[5],0)/y[5]*(yspA*(mu_maxA*y[4]/(ksA+y[4]) - kdA)/ysxA - (mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])/ysxE + msE)*y[1] - ((mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS+msS)*y[3]-saiout/V), # Ammonia
         (max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])-kdE)*y[1],0)/((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4]) - kdE)*y[1])*mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*ysp_g/ysxE + (1-max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])-kdE)*y[1],0)/(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])-kdE)*y[1])*ysp_m*msE)*y[1]] # Product

    prob=ODEProblem(f,[E,A,S,N,C,P],(0.0,tspan))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    # plot(a,A[1,:],title="Microbial concentration profile",xaxis="Time(hr)",yaxis="X(g/L)",label=false)
    # plot(a,A[2,:],title="Substrate concentration profile",xaxis="Time(hr)",yaxis="S(g/L)",label=false)
    # plot(a,A[3,:],title="Product concentration profile",xaxis="Time(hr)",yaxis="P(g/L)",label=false)
    return a,A[1,:],A[2,:],A[3,:],A[6,:],A[4,:],A[5,:]
    Et,At,St,Pt,Ct,Nt
end
