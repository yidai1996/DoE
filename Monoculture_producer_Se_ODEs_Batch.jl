# Simulation for Carbon Fixer
using Plots, JuMP, Ipopt, DifferentialEquations, NLsolve

function loadProcessData()
    global mu_max=1.7 #h^-1 from David's thesis(meeting slides from Prof.Lin) 1.7
    global ms=1 # gsubstrate/gbiomass/h +-0.0008 h^-1 from David's thesis  substrate used for maintenence
    global T0=303 #K
    global pH=7 # from Sofia medium
    global kd=0.003 #h^-1 cell death rate from David's thesis
    global ks=0.1 # gbiomass/L +-0.004 from David's thesis
    # global ki=0.3 # my guessing
    global n=4.86
    global P_star=0.20
    global ysp_g=0.3 # gbiomass/gsubstrate from Minty 13 0.322
    global ysp_m=0.2 # gbiomass/gsubstrate from David's thesis, 0.409 at the first try
    # global ysx=0.06 #http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=3 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    global ysx=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global D0= 0.68 # h^-1 initial dilusion rate
    global C0=56.4 # g/l initial substrate feeding concentration
    global X0=0.01 # g/L initial cell concentration 0.7g/L is from Figure 2.4 on Page 34 of David's thesis
    global S0=1 # g/L initial substrate concentration
    global P0=0.00 # g/L initial product concentration
    # global Xs=0.7 # g/L set point of cell concentration
    # global Ps=0.25 # g/L
    # global Ss=5.0
    global tspan=100
    println("Parameters Loaded!")
end

function EcoliGrowth(X0,S0,P0,tspan)
    global tt1,X1,S1,P1=ODEStep(X0,S0,P0,tspan)
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
    # plot(tt1,X1/0.396,title="Microbial concentration profile",xaxis="Time(hr)",yaxis="OD600",label=false)
    # savefig("0316_X_soln_S0_2.png")
    # plot(tt1,dxdt,title="Microbial growth rate",xaxis="Time(hr)",yaxis="dXdt(g/L/h)",label=false)
    # plot(tt1,S1,title="Substrate concentration profile",xaxis="Time(hr)",yaxis="S(g/L)",label=false)
    # savefig("0316_S_soln_batch_P_star_0.05.png")
    # plot(tt1,dSdt,title="Substrate change rate",xaxis="Time(hr)",yaxis="dSdt(g/L/h)",label=false)
    # plot(tt1,P1,title="Product concentration profile",xaxis="Time(hr)",yaxis="P(g/L)",label=false)
    # savefig("0316_P_soln_batch_P_star_0.05.png")
    # plot(tt1,dPdt,title="Product rate",xaxis="Time(hr)",yaxis="dPdt(g/L/h)",label=false)
end

function ODEStep(X,S,P,tspan) # Use one ODE solver to solve the whole system
    f(y,p,t)=[(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2]) - kd)*y[1],
         -max(y[2],0)/y[2]*(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])/ysx + ms)*y[1],
         0]
         # (max((mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])-kd)*y[1],0)/((mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2]) - kd)*y[1])*mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*ysp_g/ysx + (1-max((mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])-kd)*y[1],0)/(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])-kd)*y[1])*ysp_m*ms)*y[1]] # X,S,P
    prob=ODEProblem(f,[X,S,P],(0.0,tspan))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    # plot(a,A[1,:],title="Microbial concentration profile",xaxis="Time(hr)",yaxis="X(g/L)",label=false)
    # plot(a,A[2,:],title="Substrate concentration profile",xaxis="Time(hr)",yaxis="S(g/L)",label=false)
    # plot(a,A[3,:],title="Product concentration profile",xaxis="Time(hr)",yaxis="P(g/L)",label=false)
    return a,A[1,:],A[2,:],A[3,:]
end
