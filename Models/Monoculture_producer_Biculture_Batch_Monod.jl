# Simulation for Continuously Triculturing System
# Simplest version: the limited substrates are sucrose and ammonia for Av/Se respectively, while E.coli consider sucrose and oxygen transfer limitation
using Plots, DifferentialEquations, NLsolve, XLSX, Printf

function loadProcessData()
    # Fitted parameters:
    global mu_maxA=0.11
    global mu_maxS=0.0217
    global msA=0.05
    global msS=0.05
    global ksA=0.520
    global ksS=0.344
    global ysxA=2.451
    global ysxS=2.54
    global kdA=0.01
    global kdS=0.1
    global yspA=0.018
    global mu_maxE=0.34 #h^-1 from David's thesis(meeting slides from Prof.Lin) 1.7
    global msE=0.1 # gsubstrate/gbiomass/h +-0.0008 h^-1 from David's thesis  substrate used for maintenence
    global msE_NH4=0.1 # gsubstrate/gbiomass/h
    global T0=303 #K
    global pH=7 # from Sofia medium
    global kdE=0.003 #h^-1 cell death rate from David's thesis
    global ksE=0.1 # gbiomass/L +-0.004 from David's thesis
    global ksE_NH4=0.1
    global n=4.86
    global P_star=0.20
    global ysp_g=0.3 # 
    global ysp_m=0.2 # 
    global ysxE=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    global ysxE_NH4=1.017
    global ysxA=0.17 # https://www.researchgate.net/figure/Growth-kinetics-of-Azotobacter-vinelandii-in-medium-before-and-after-optimization_tbl2_301753155
    global ysxS=0.17 # https://www.sciencedirect.com/science/article/pii/S0960852406004792
    global A0=1 # g/l initial substrate feeding concentration
    global E0=0.01 # g/L initial cell concentration 0.7g/L is from Figure 2.4 on Page 34 of David's thesis   
    global S0=0.5 # g/L initial substrate concentration
    # global N0=1 # g/L initial substrate concentration
    global N0=100 # g/L initial substrate concentration
    # global C0=1 # g/L initial substrate concentration
    global C0=30 # g/L initial substrate concentration
    global P0=0.001 # g/L initial product concentration
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
    global saiin=0.05  # L/min
    global saioutE=0.01 # David's thesis
    global saioutS=0.01
    global saioutA=0.01
    global V=0.5 #L
    global Vt=3*V
    global DE=saioutE/V # 0.2
    global DA=saioutA/V # 0.2
    global DS=saioutS/V # 0.2
    global D0=[DE DA DS]
    global KA=300
    global KS=300
    global out_dir="G:\\My Drive\\Research\\DOE project\\Modeling\\GitDoE"
    # println("Parameters Loaded!")
end

function CocultureGrowth(parameterS,parameterA,initialC) # Continuous flow
    println(parameterS,parameterA,initialC)
    # parameterS=[mu_maxS ksS ysxS yspS kdS]
    # parameterA=[mu_maxA ksA ysxA yspA kdA]
    # initialC = [Av0 Se0 C0 N0]
    loadProcessData()
    global tt1,At1,St1,Ct1,Nt1=Startup(initialC[1],initialC[2],initialC[3],initialC[4],tspan1,parameterS,parameterA)
    # Startup stage
    # global dA1dt=zeros(size(tt1)[1])
    # global dS1dt=zeros(size(tt1)[1])
    # for i=1:size(tt1)[1]
    #     dA1dt[i]=(mu_maxA*Ct1[i]/(ksA+Ct1[i]) - kdA)*At1[i]
    #     # dA1dt[i]=mu_maxA*(KA-At1[i])/KA*At1[i]
    # end
    # # println(size(tt1)[1])
    # for i=1:size(tt1)[1]
    #     dS1dt[i]=(mu_maxS*Nt1[i]/(ksS+Nt1[i]) - kdS)*St1[i]
    #     # dS1dt[i]=mu_maxS*(KS-St1[i])/KS*St1[i]
    # end
    s1=parameterS[1]
    s2=parameterS[2]
    s3=parameterS[3]
    s4=parameterS[4]
    s5=parameterS[5]
    a1=parameterA[1]
    a2=parameterA[2]
    a3=parameterA[3]
    a4=parameterA[4]
    a5=parameterA[5]
    c1=initialC[1]
    c2=initialC[2]
    c3=initialC[3]
    c4=initialC[4]
    # plot(tt1,dA1dt,label="Growth rate of Av when biculture",xaxis="Time(hr)",yaxis="g/L/h",title="Start up profiles for co-culture",framestyle=:box,legend=:topleft)
    # plot!(tt1,dS1dt,label="Growth rate of Se when biculture",xaxis="Time(hr)",yaxis="g/L/h",)
    # savefig("muS_$(@sprintf("%.2f",s1)) ksS_$(@sprintf("%.2f",s2)) ysxS_$(@sprintf("%.2f",s3)) yspS_$(@sprintf("%.2f",s4)) msS_$(@sprintf("%.2f",s5)) muA_$(@sprintf("%.2f",a1)) ksA_$(@sprintf("%.2f",a2)) ysxA_$(@sprintf("%.2f",a3)) yspA_$(@sprintf("%.2f",a4)) msA_$(@sprintf("%.2f",a5)).pdf")
 
    # Microbial
    plot(tt1,At1,label="Av concentration profile when biculture",xaxis="Time(hr)",yaxis="Av(g/L)",framestyle=:box,legend=:topleft)
    plot!(tt1,St1,label="Se concentration profile when biculture",xaxis="Time(hr)",yaxis="Se(g/L)")
    savefig("Microbial_muS_$(@sprintf("%.2f",s1)) ksS_$(@sprintf("%.2f",s2)) ysxS_$(@sprintf("%.2f",s3)) yspS_$(@sprintf("%.2f",s4)) msS_$(@sprintf("%.2f",s5)) muA_$(@sprintf("%.2f",a1)) ksA_$(@sprintf("%.2f",a2)) ysxA_$(@sprintf("%.2f",a3)) yspA_$(@sprintf("%.2f",a4)) msA_$(@sprintf("%.2f",a5)) A0_$(@sprintf("%.2f",c1)) S0_$(@sprintf("%.2f",c2)).pdf")
    plot(tt1,Ct1,label="Sucrose concentration profile when biculture",xaxis="Time(hr)",yaxis="Sucrose(g/L)",framestyle=:box,legend=:topleft)
    plot!(tt1,Nt1,label="Ammonia concentration profile when biculture",xaxis="Time(hr)",yaxis="Ammonia(g/L)")
    savefig("Nutrient_muS_$(@sprintf("%.2f",s1)) ksS_$(@sprintf("%.2f",s2)) ysxS_$(@sprintf("%.2f",s3)) yspS_$(@sprintf("%.2f",s4)) msS_$(@sprintf("%.2f",s5)) muA_$(@sprintf("%.2f",a1)) ksA_$(@sprintf("%.2f",a2)) ysxA_$(@sprintf("%.2f",a3)) yspA_$(@sprintf("%.2f",a4)) msA_$(@sprintf("%.2f",a5)).pdf")

    # Store data into excel files
    println("writing plots to files")
    println(out_dir)
    top_excel_file = out_dir * "\\muS_$(@sprintf("%.2f",s1)) ksS_$(@sprintf("%.2f",s2)) ysxS_$(@sprintf("%.2f",s3)) yspS_$(@sprintf("%.2f",s4)) msS_$(@sprintf("%.2f",s5)) muA_$(@sprintf("%.2f",a1)) ksA_$(@sprintf("%.2f",a2)) ysxA_$(@sprintf("%.2f",a3)) yspA_$(@sprintf("%.2f",a4)) msA_$(@sprintf("%.2f",a5)).xlsx"
    column_names = ["times (hr)","Av","Se", "Sucrose", "Ammonia"]
    data=[tt1,At1,St1,Ct1,Nt1]
    # write to excel file
    XLSX.writetable(top_excel_file, data, column_names)
  
end

function Startup(A,S,N,C,tspan1,pS,pA) # batch
    # pS=[mu_maxS ksS ysxS yspS msS]
    # pA=[mu_maxA ksA ysxA yspA msA]
    # Monod Equation
    f(y,p,t)=[(pA[1]*y[3]/(pA[2]+y[3]) - kdA)*y[1],# X(Av)
         (pS[1]*y[4]/(pS[2]+y[4]) - kdS)*y[2],# X(Se)
         max(y[3],0)/y[3]*( - (pA[1]/(pA[2]+y[3])/pA[3]+pA[5])*y[1] + pS[4]*(pS[1]*y[4]/(ksS+y[4]))/pS[3]*y[2]), # Sucrose
         max(y[4],0)/y[4]*(pA[4]*pA[1]*y[3]/(pA[2]+y[3])*y[1]/pA[3]  - (pS[1]*y[4]/(pS[2]+y[4])/pS[3]+pS[5])*y[2])] # Ammonia
    prob=ODEProblem(f,[A,S,N,C],(0.0,tspan1))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    return a,A[1,:],A[2,:],A[3,:],A[4,:]
    # At,St,Ct,Nt
end

function StartupNoms(A,S,N,C,tspan1,pS,pA) # batch
    # pS=[mu_maxS ksS ysxS yspS kdS]
    # pA=[mu_maxA ksA ysxA yspA kdA]
    # Monod Equation
    f(y,p,t)=[(pA[1]*y[3]/(pA[2]+y[3]) - pA[5])*y[1],# X(Av)
         (pS[1]*y[4]/(pS[2]+y[4]) - pS[5])*y[2],# X(Se)
         max(y[3],0)/y[3]*( - (pA[1]/(pA[2]+y[3])/pA[3]+pA[5])*y[1] + pS[4]*(pS[1]*y[4]/(ksS+y[4]))/pS[3]*y[2]), # Sucrose
         max(y[4],0)/y[4]*(pA[4]*pA[1]*y[3]/(pA[2]+y[3])*y[1]/pA[3]  - (pS[1]*y[4]/(pS[2]+y[4])/pS[3]+pS[5])*y[2])] # Ammonia
    prob=ODEProblem(f,[A,S,N,C],(0.0,tspan1))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    return a,A[1,:],A[2,:],A[3,:],A[4,:]
    # At,St,Ct,Nt
end

function Tripartite(D,E,A,S,N,C,P,tspan2) # Use one ODE solver to solve the whole system
    # D: dilusion rate
    # E/A/S: e.coli/Ntrigen/Carbon fixer  concentration
    # With inhibition from isobutanol
    # f(y,p,t)=[(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1]-D[1]*y[1], # X(E.coli)
    #      (mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]-D[2]*y[2],# X(Av)
    #      (mu_maxS*y[5]/(ksS+y[5]) - kdS)*y[3]-D[3]*y[3],# X(Se)
    #      1/Vt*max(y[4],0)/y[4]*(-(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE + msE)*y[1] - ((mu_maxA*y[4]/(ksA+y[4]) - kdA)/ysxA+msA)*y[2] + yspS*(mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS*y[3]  - sum(D)*y[4]), # Sucrose
    #      1/Vt*max(y[5],0)/y[5]*(yspA*(mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]/ysxA - (mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE_NH4 + msE_NH4)*y[1] - ((mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS+msS)*y[3] - sum(D)*y[5]), # Ammonia
    #      1/Vt*(max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1])*mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])*ysp_g/ysxE + (1-max((mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/(mu_maxE*y[4]*(max(1-y[6]/P_star,0.0))^n/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1])*ysp_m*msE)*y[1] - sum(D)*y[6]] # Product

    # Without inhibition
    f(y,p,t)=[(mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1]-D[1]*y[1], # X(E.coli)
         (mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]-D[2]*y[2],# X(Av)
         (mu_maxS*y[5]/(ksS+y[5]) - kdS)*y[3]-D[3]*y[3],# X(Se)
         1/Vt*max(y[4],0)/y[4]*(-(mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE + msE)*y[1] - ((mu_maxA*y[4]/(ksA+y[4]) - kdA)/ysxA+msA)*y[2] + yspS*(mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS*y[3]  - sum(D)*y[4]), # Sucrose
         1/Vt*max(y[5],0)/y[5]*(yspA*(mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]/ysxA - (mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE_NH4 + msE_NH4)*y[1] - ((mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS+msS)*y[3] - sum(D)*y[5]), # Ammonia
         1/Vt*(max((mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/((mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1])*mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])*ysp_g/ysxE + (1-max((mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/((mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1]))*ysp_m*msE)*y[1] - sum(D)*y[6]] # Product

    prob=ODEProblem(f,[E,A,S,N,C,P],(0.0,tspan2))
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
