# Simulation for Continuously Triculturing System
# Simplest version: the limited substrates are sucrose and ammonia for Av/Se respectively, while E.coli consider sucrose and oxygen transfer limitation
using Plots, JuMP, Ipopt, DifferentialEquations, NLsolve, XLSX

function loadProcessData()
    # Fitted parameters:
    global mu_maxA=0.11
    global mu_maxS=0.0217
    global msA=4
    global msS=0.79
    global ksA=0.520
    global ksS=0.344
    global ysxA=2.451
    global ysxS=2.54
    global kdA=0.006
    global kdS=0
    global yspA=0.018
    global mu_maxE=0.34 #h^-1 from David's thesis(meeting slides from Prof.Lin) 1.7
    # global mu_maxA=0.34 #h^-1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5478974/
    # global mu_maxA=0.091 #h^-1 https://www.sciencedirect.com/science/article/pii/S1369703X02001766
    # global mu_maxS=0.0504 #h^-1 https://onlinelibrary.wiley.com/doi/full/10.1002/cjce.22154
    # global mu_maxS=0.34 #h^-1 https://onlinelibrary.wiley.com/doi/full/10.1002/cjce.22154
    global msE=0.1 # gsubstrate/gbiomass/h +-0.0008 h^-1 from David's thesis  substrate used for maintenence
    global msE_NH4=0.1 # gsubstrate/gbiomass/h
    # global msA=0.31 # gsubstrate/gbiomass/h
    # global msA=0.22 # gsubstrate/gbiomass/h
    # global msS=0.22 # gsubstrate/gbiomass/h
    global T0=303 #K
    global pH=7 # from Sofia medium
    global kdE=0.003 #h^-1 cell death rate from David's thesis
    # global kdA=0.003 #h^-1
    # global kdS=0.003 #h^-1
    global ksE=0.1 # gbiomass/L +-0.004 from David's thesis
    # global ksA=0.1 # gbiomass/L https://www.sciencedirect.com/science/article/pii/S1369703X02001766
    # global ksS=0.1 # gbiomass/L https://onlinelibrary.wiley.com/doi/full/10.1111/j.1529-8817.2005.04063.x
    global ksE_NH4=0.1
    # global ki=0.3 # my guessing
    global n=4.86
    global P_star=0.20
    global ysp_g=0.3 # 
    global ysp_m=0.2 # 
    # global yspA=2 # https://microbialcellfactories.biomedcentral.com/track/pdf/10.1186/s12934-020-01362-9.pdf
    # global yspS=2 #
    # global ysx=0.06 #http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    # global ysx=3 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    global ysxE=1.017 # http://staff.du.edu.eg/upfilestaff/1066/researches/31066_1619277717__jawed2020._.pdf
    global ysxE_NH4=1.017
    global ysxA=0.17 # https://www.researchgate.net/figure/Growth-kinetics-of-Azotobacter-vinelandii-in-medium-before-and-after-optimization_tbl2_301753155
    global ysxS=0.17 # https://www.sciencedirect.com/science/article/pii/S0960852406004792
    # global D0= 0.68 # h^-1 initial dilusion rate
    global A0=0.01 # g/l initial substrate feeding concentration
    global E0=0.01 # g/L initial cell concentration 0.7g/L is from Figure 2.4 on Page 34 of David's thesis
    global S0=0.01 # g/L initial substrate concentration
    global N0=1 # g/L initial substrate concentration
    global C0=1 # g/L initial substrate concentration
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
    global tspan1=48 # h David's thesis
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
    global out_dir="G:\\My Drive\\Research\\DOE project\\Modeling\\Triculture\\modified triculture model\\triculture test without inhibition from isobutanol (Monod equation)"
    println("Parameters Loaded!")
end

function CocultureGrowth() # Continuous flow
    loadProcessData()
    global tt1,At1,St1,Ct1,Nt1=Startup(A0,S0,N0,C0,tspan1)
    # Startup stage
    global dA1dt=zeros(size(tt1)[1])
    global dS1dt=zeros(size(tt1)[1])
    for i=1:size(tt1)[1]
        dA1dt[i]=(mu_maxA*Ct1[i]/(ksA+Ct1[i]) - kdA)*At1[i]
        # dA1dt[i]=mu_maxA*(KA-At1[i])/KA*At1[i]
    end
    println(size(tt1)[1])
    for i=1:size(tt1)[1]
        dS1dt[i]=(mu_maxS*Nt1[i]/(ksS+Nt1[i]) - kdS)*St1[i]
        # dS1dt[i]=mu_maxS*(KS-St1[i])/KS*St1[i]
    end
    plot(tt1,dA1dt,label="Growth rate of Av when biculture",xaxis="Time(hr)",yaxis="g/L/h",title="Start up profiles for co-culture",framestyle=:box,legend=:topleft)
    plot!(tt1,dS1dt,label="Growth rate of Se when biculture",xaxis="Time(hr)",yaxis="g/L/h",)
    savefig("Coculture profile monod equation of growth rate tspan1_48.pdf")
    plot(tt1,At1,label="Av concentration profile when biculture",xaxis="Time(hr)",yaxis="Av(g/L)",framestyle=:box,legend=:topleft)
    plot!(tt1,St1,label="Se concentration profile when biculture",xaxis="Time(hr)",yaxis="Se(g/L)")
    savefig("Coculture profile monod equation of microbial profile tspan1_48.pdf")
    plot(tt1,Ct1,label="Sucrose concentration profile when biculture",xaxis="Time(hr)",yaxis="Sucrose(g/L)",framestyle=:box,legend=:topleft)
    plot!(tt1,Nt1,label="Ammonia concentration profile when biculture",xaxis="Time(hr)",yaxis="Ammonia(g/L)")
    savefig("Coculture profile monod equation of nutrient tspan1_48.pdf")

    # Store data into excel files
    # println("writing plots to files")
    # top_excel_file = out_dir * "\\Profiles of All Microbial without inhibition tspan_48.xlsx"
    # column_names = ["times (hr)","Av","Se", "Sucrose", "Ammonia","Growth rate of Av","Growth rate of Se"]
    # data=[tt1,At1,St1,Ct1,Nt1,dA1dt,dS1dt]
    # # write to excel file
    # XLSX.writetable(top_excel_file, data, column_names)
  
end

function AllGrowth() # Continuous flow
    loadProcessData()
    global tt1,At1,St1,Ct1,Nt1=Startup(A0,S0,N0,C0,tspan1)
    # Startup stage
    global dA1dt=zeros(size(tt1)[1])
    global dS1dt=zeros(size(tt1)[1])
    for i=1:size(tt1)[1]
        dA1dt[i]=(mu_maxA*Ct1[i]/(ksA+Ct1[i]) - kdA)*At1[i]
    end
    println(size(tt1)[1])
    for i=1:size(tt1)[1]
        dS1dt[i]=(mu_maxS*Nt1[i]/(ksS+Nt1[i]) - kdS)*St1[i]
    end
    # plot(tt1,dA1dt,label="Growth rate of Av when biculture",xaxis="Time(hr)",yaxis="g/L/h",title="Start up profiles for co-culture",framestyle=:box,legend=:topleft)
    # plot!(tt1,dS1dt,label="Growth rate of Se when biculture",xaxis="Time(hr)",yaxis="g/L/h",)
    # plot!(tt1,At1,label="Av concentration profile when biculture",xaxis="Time(hr)",yaxis="Av(g/L)")
    # plot!(tt1,St1,label="Se concentration profile when biculture",xaxis="Time(hr)",yaxis="Se(g/L)")
    # plot!(tt1,Ct1,label="Sucrose concentration profile when biculture",xaxis="Time(hr)",yaxis="Sucrose(g/L)")
    # plot!(tt1,Nt1,label="Ammonia concentration profile when biculture",xaxis="Time(hr)",yaxis="Ammonia(g/L)")
    # savefig("Coculture profile (with ms parameters) tspan1_36.pdf")

    # Incubate E.coli to start triculturing
    global tt,Et,At,St,Pt,Ct,Nt=Tripartite(D0,E0,At1[end],St1[end],Ct1[end],Nt1[end],P0,tspan2)
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
    size_t=size(tt)[1]+size(tt1)[1]
    global final_dEdt=zeros(size_t)
    global final_dAdt=zeros(size_t)
    global final_dSdt=zeros(size_t)
    global final_Et=zeros(size_t)
    global final_At=zeros(size_t)
    global final_St=zeros(size_t)
    global final_Ct=zeros(size_t)
    global final_Nt=zeros(size_t)
    global final_Pt=zeros(size_t)
    for i=1:size(tt)[1]
        dEdt[i]=(mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i]) - kdE)*Et[i]-saioutE*Et[i]
    end
    
    ttt=cat(tt1,tt.+tspan1;dims=(1,1))
    println(size(ttt)[1])
    # dE1dt=zeros(size(tt1)[1])
    # final_dEdt=cat(dE1dt,dEdt;dims=(1,1))
    # println(size(ttt))
    # println(size(final_dEdt))
    # plot(ttt,final_dEdt,title="Growth rate of E.coli",xaxis="Time(hr)",yaxis="g/L/h",label=false)
    # savefig("E.coliGrowthRate_test.png")

    # for i=1:size(tt)[1]
    #     dAdt[i]=(mu_maxA*Ct[i]/(ksA+Ct[i]) - kdA)*At[i] - saioutA*At[i]
    # end
    # final_dAdt=cat(dA1dt,dAdt;dims=(1,1))
    # plot(ttt,final_dAdt,title="Growth rate of Av",xaxis="Time(hr)",yaxis="g/L/h",label=false)
    # savefig("AvGrowthRate_test.png")

    # for i=1:size(tt)[1]
    #     dSdt[i]=(mu_maxS*Nt[i]/(ksS+Nt[i]) - kdS)*St[i] - saioutS*St[i]
    # end
    # final_dSdt=cat(dS1dt,dSdt;dims=(1,1))
    # plot(ttt,final_dSdt,title="Growth rate of Se",xaxis="Time(hr)",yaxis="g/L/h",label=false)
    # savefig("SeGrowthRate_test.png")

    # for i=1:size(tt)[1]
    #     dCdt1[i]= - (mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i])/ysxE + msE)*Et[i]
    #     dCdt2[i]= - ((mu_maxA*Ct[i]/(ksA+Ct[i]) - kdA)/ysxA+msA)*At[i]
    #     dCdt3[i]= + yspS*(mu_maxS*Nt[i]/(ksS+Nt[i]) - kdS)/ysxS
    # end
    # plot(tt,dCdt1,title="Sucrose consumed/produced rate",xaxis="Time(hr)",yaxis="g/L/h",label="Consumed by E.coli")
    # plot!(tt,dCdt2,xaxis="Time(hr)",yaxis="g/L/h",label="Consumed by Av")
    # plot!(tt,dCdt3,xaxis="Time(hr)",yaxis="g/L/h",label="Produced by Se")
    # savefig("Sucrose consumedandproduced rate_test.png")

    # for i=1:size(tt)[1]
    #     dNdt1[i]= yspA*(mu_maxA*Ct[i]/(ksA+Ct[i]) - kdA)*At[i]/ysxA
    #     dNdt2[i]= - (mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i])/ysxE + msE)*Et[i]
    #     dNdt3[i]= - ((mu_maxS*Nt[i]/(ksS+Nt[i]) - kdS)/ysxS+msS)*St[i]
    # end
    # plot(tt,dNdt1,title="Ammonia consumed/produced rate",xaxis="Time(hr)",yaxis="g/L/h",label="Produced by Av")
    # plot!(tt,dNdt2,xaxis="Time(hr)",yaxis="g/L/h",label="Consumed by E.coli")
    # plot!(tt,dNdt3,xaxis="Time(hr)",yaxis="g/L/h",label="Consumed by Se")
    # savefig("Ammonia consumedandproduced rate_test.png")

    # for i=1:size(tt)[1]
    #     zz=(mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i]) - kdE)*Et[i]
    #     dPdt[i]=(max(zz,0)/zz*mu_maxE*Ct[i]*(max(1-Pt[i]/P_star,0.0))^n/(ksE+Ct[i])*Nt[i]/(ksE_NH4+Nt[i])*ysp_g/ysxE + (1-max(zz,0)/zz)*ysp_m*msE)*Et[i]
    # end
    # plot(tt,dPdt,title="Isobutanol produced rate",xaxis="Time(hr)",yaxis="g/L/h",label=false)
    # savefig("Isobutanol produced rate_test.png")

    final_Et=cat(zeros(size(tt1)[1]),Et;dims=(1,1))
    plot(ttt,final_Et/0.396,title="E.coli concentration profile",xaxis="Time(hr)",yaxis="OD600",label=false)
    savefig("Ecoli_test triculture.png")

    final_At=cat(At1,At;dims=(1,1))
    plot(ttt,final_At,title="Av concentration profile",xaxis="Time(hr)",yaxis="Av(g/L)",label=false)
    savefig("Av_test triculture.png")

    final_St=cat(St1,St;dims=(1,1))
    plot(ttt,final_St,title="Se concentration profile",xaxis="Time(hr)",yaxis="Se(g/L)",label=false)
    savefig("Se_test triculture.png")

    final_Ct=cat(Ct1,Ct;dims=(1,1))
    plot(ttt,final_Ct,title="Sucrose concentration profile",xaxis="Time(hr)",yaxis="Sucrose(g/L)",label=false)
    savefig("Sucrose_test triculture.png")

    final_Nt=cat(Nt1,Nt;dims=(1,1))
    plot(ttt,final_Nt,title="Ammonia concentration profile",xaxis="Time(hr)",yaxis="Ammonia(g/L)",label=false)
    savefig("Ammonia_test triculture.png")

    final_Pt=cat(zeros(size(tt1)[1]),Pt;dims=(1,1))
    plot(ttt,final_Pt,title="Product concentration profile",xaxis="Time(hr)",yaxis="Isobutanol(g/L)",label=false)
    savefig("Isobutanol_test.png")

    # Store data into excel files
    println("writing plots to files")
    top_excel_file = out_dir * "\\Profiles of All Microbial without inhibition.xlsx"
    column_names = ["times (hr)","Ec","Av","Se", "Sucrose", "Ammonia","Isobutanol"]
    data=[ttt,final_Et,final_At,final_St,final_Ct,final_Nt,final_Pt]
    # write to excel file
    XLSX.writetable(top_excel_file, data, column_names)
  
end

function Startup(A,S,N,C,tspan1) # batch
    # Logistic equation
    # f(y,p,t)=[mu_maxA*(KA-y[1])/KA*y[1]-kdA*y[1],# X(Av) with death rate
    # mu_maxS*(KS-y[2])/KS*y[2]-kdS*y[2],# X(Se)
    # # f(y,p,t)=[mu_maxA*(KA-y[1])/KA*y[1],# X(Av)
    # #      mu_maxS*(KS-y[2])/KS*y[2],# X(Se)
    #      max(y[3],0)/y[3]*( - (mu_maxA*(KA-y[1])/KA*y[1]/ysxA+msA)*y[1] + yspS*mu_maxS*(KS-y[2])/KS*y[2]/ysxS*y[2]), # Sucrose
    #      max(y[4],0)/y[4]*(yspA*mu_maxA*(KA-y[1])/KA*y[1]/ysxA  - (mu_maxS*(KS-y[2])/KS*y[2]/ysxS+msS)*y[2])] # Ammonia
    # Monod Equation
    f(y,p,t)=[(mu_maxA*y[3]/(ksA+y[3]) - kdA)*y[1],# X(Av)
         (mu_maxS*y[4]/(ksS+y[4]) - kdS)*y[2],# X(Se)
         max(y[3],0)/y[3]*( - (mu_maxA*y[3]/(ksA+y[3])/ysxA)*y[1] + yspS*(mu_maxS*y[4]/(ksS+y[4]))/ysxS*y[2]), # Sucrose
         max(y[4],0)/y[4]*(yspA*mu_maxA*y[3]/(ksA+y[3])*y[1]/ysxA  - (mu_maxS*y[4]/(ksS+y[4])/ysxS)*y[2])] # Ammonia
        #  max(y[3],0)/y[3]*( - (mu_maxA*y[3]/(ksA+y[3])/ysxA+msA)*y[1] + yspS*(mu_maxS*y[4]/(ksS+y[4]))/ysxS*y[2]), # Sucrose
        #  max(y[4],0)/y[4]*(yspA*mu_maxA*y[3]/(ksA+y[3])*y[1]/ysxA  - (mu_maxS*y[4]/(ksS+y[4])/ysxS+msS)*y[2])] # Ammonia
    prob=ODEProblem(f,[A,S,N,C],(0.0,tspan1))
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
