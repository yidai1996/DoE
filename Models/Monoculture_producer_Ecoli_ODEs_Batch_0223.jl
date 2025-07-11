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
    global tspan1=72 # h David's thesis
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

# function loadProcessData()
#     global mu_max=1.7 #h^-1 from David's thesis(meeting slides from Prof.Lin) 1.7
#     global ms=1 # gsubstrate/gbiomass/h +-0.0008 h^-1 from David's thesis  substrate used for maintenence
#     global T0=303 #K
#     global pH=7 # from Sofia medium
#     global kd=0.003 #h^-1 cell death rate from David's thesis
#     global ks=0.1 # gbiomass/L +-0.004 from David's thesis
#     # global ki=0.3 # my guessing
#     global n=4.86
#     global P_star=0.20
#     global ysp_g=0.3 # gbiomass/gsubstrate from Minty 13 0.322
#     global ysp_m=0.2 # gbiomass/gsubstrate from David's thesis, 0.409 at the first try
#     # global ysx=0.06 #http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
#     # global ysx=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
#     # global ysx=3 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
#     global ysx=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
#     # global D0= 0.68 # h^-1 initial dilusion rate
#     global C0=56.4 # g/l initial substrate feeding concentration
#     global X0=0.01 # g/L initial cell concentration 0.7g/L is from Figure 2.4 on Page 34 of David's thesis
#     global S0=1 # g/L initial substrate concentration
#     global P0=0.00 # g/L initial product concentration
#     # global Xs=0.7 # g/L set point of cell concentration
#     # global Ps=0.25 # g/L
#     # global Ss=5.0
#     global kO2=0.1 # guessing
#     global M_O2=16 # g/mol
#     global yox=1.017
#     global mo=0.01 # guessing
#     global kLaO=2.766*60 # 1/h # https://www.sciencedirect.com/science/article/pii/S0032959200002727
#     global xO2_sat=7.5/16/1000 # mol/L https://www.waterboards.ca.gov/water_issues/programs/swamp/docs/cwt/guidance/3110en.pdf
#     global O20=xO2_sat
#     global tspan=100 # h
#     println("Parameters Loaded!")
# end

function EcoliGrowth()
    loadProcessData()
    # global tt1,X1,S1,P1=ODEStep(X0,S0,P0,tspan)
    global tt,Xt,St,Pt,O2=ODEStep(X0,S0,P0,O20,tspan)
    plot(tt1,X1/0.396,title="Microbial concentration profile",xaxis="Time(hr)",yaxis="OD600",label=false)
    # savefig("0316_X_soln_S0_2.png")
    # plot(tt1,S1,title="Substrate concentration profile",xaxis="Time(hr)",yaxis="S(g/L)",label=false)
    # savefig("0316_S_soln_batch_P_star_0.05.png")
    # plot(tt1,P1,title="Product concentration profile",xaxis="Time(hr)",yaxis="P(g/L)",label=false)
    # savefig("0316_P_soln_batch_P_star_0.05.png")

    # Calculate growth/consuming rate
    # global dxdt=zeros(size(tt1)[1])
    # global dPdt=zeros(size(tt1)[1])
    # global dSdt=zeros(size(tt1)[1])
    # global dOdt=zeros(size(tt1)[1])
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
    # plot(tt1,dxdt,title="Microbial growth rate",xaxis="Time(hr)",yaxis="dXdt(g/L/h)",label=false)
    # plot(tt1,dSdt,title="Substrate change rate",xaxis="Time(hr)",yaxis="dSdt(g/L/h)",label=false)
    # plot(tt1,dPdt,title="Product rate",xaxis="Time(hr)",yaxis="dPdt(g/L/h)",label=false)
end

function ODEStep(X,S,P,O2,tspan) # Use one ODE solver to solve the whole system
    # f(y,p,t)=[(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2]) - kd)*y[1],
    #      -max(y[2],0)/y[2]*(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])/ysx + ms)*y[1],
    #      0] # No production
    # f(y,p,t)=[(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2]) - kd)*y[1],
    #           -max(y[2],0)/y[2]*(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])/ysx + ms)*y[1],
    #          (max((mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])-kd)*y[1],0)/((mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2]) - kd)*y[1])*mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*ysp_g/ysx + (1-max((mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])-kd)*y[1],0)/(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])-kd)*y[1])*ysp_m*ms)*y[1]]
    #          # Producing isobutanol
    f(y,p,t)=[(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*M_O2*y[4]/(kO2+M_O2*y[4]) - kd)*y[1],
               -max(y[2],0)/y[2]*(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*M_O2*y[4]/(kO2+M_O2*y[4])/ysx + ms)*y[1],
              (max((mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*M_O2*y[4]/(kO2+M_O2*y[4])-kd)*y[1],0)/((mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*M_O2*y[4]/(kO2+M_O2*y[4]) - kd)*y[1])*mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*M_O2*y[4]/(kO2+M_O2*y[4])*ysp_g/ysx + (1-max((mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*M_O2*y[4]/(kO2+M_O2*y[4])-kd)*y[1],0)/(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*M_O2*y[4]/(kO2+M_O2*y[4])-kd)*y[1])*ysp_m*ms)*y[1],
              -max(y[4],0)/y[4]*(mu_max*y[2]*(max(1-y[3]/P_star,0.0))^n/(ks+y[2])*M_O2*y[4]/(kO2+M_O2*y[4])/yox+mo)*y[1] + M_O2*kLaO*(xO2_sat-y[4])]
              # Considering oxygen transfer limitation in Monod equation with production
    prob=ODEProblem(f,[X,S,P,O2],(0.0,tspan))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    # plot(a,A[1,:],title="Microbial concentration profile",xaxis="Time(hr)",yaxis="X(g/L)",label=false)
    # plot(a,A[2,:],title="Substrate concentration profile",xaxis="Time(hr)",yaxis="S(g/L)",label=false)
    # plot(a,A[3,:],title="Product concentration profile",xaxis="Time(hr)",yaxis="P(g/L)",label=false)
    return a,A[1,:],A[2,:],A[3,:],A[4,:]
end
