# Simulations of Dimensionless ODEs
using Plots, DifferentialEquations, NLsolve, XLSX, Printf

function loadProcessData()
    global mu_maxA=0.1
    global mu_maxS=0.1
    # global mu_maxA=0.11
    # global mu_maxS=0.0217
    global kdA=0.025
    global kdS=0.005
    global N1 = 1
    tc = N1/mu_maxS
    global N2 = tc*mu_maxA
    global N3 = tc*kdS
    global N4 = tc*kdA
    global N5 = 3.5
    global N6 = 1
    global N7 = 1
    global N8 = 1
    global N9 = 3
    global N10 = 1

    global A0 = 0.01
    global S0 = 0.01
    global tspan1=800
    global out_dir="G:\\My Drive\\Research\\DOE project\\Modeling\\DimensionlessAnalysis\\"

    println("Parameters Loaded!")
    return [N1 N3 N5 N7 N9], [N2 N4 N6 N8 N10]
end

function CocultureGrowth(N0, C0; filename = "NotSpecific") # Continuous flow
    parameterS,parameterA = loadProcessData()
    global tt1,At1,St1,Ct1,Nt1=Startup(A0,S0,N0,C0,tspan1,parameterS,parameterA)
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
    # plot(tt1,dA1dt,label="Growth rate of Av when biculture",xaxis="Time(hr)",yaxis="g/L/h",title="Start up profiles for co-culture",framestyle=:box,legend=:topleft)
    # plot!(tt1,dS1dt,label="Growth rate of Se when biculture",xaxis="Time(hr)",yaxis="g/L/h",)
    # savefig("muS_$(@sprintf("%.2f",s1)) ksS_$(@sprintf("%.2f",s2)) ysxS_$(@sprintf("%.2f",s3)) yspS_$(@sprintf("%.2f",s4)) msS_$(@sprintf("%.2f",s5)) muA_$(@sprintf("%.2f",a1)) ksA_$(@sprintf("%.2f",a2)) ysxA_$(@sprintf("%.2f",a3)) yspA_$(@sprintf("%.2f",a4)) msA_$(@sprintf("%.2f",a5)).pdf")
 
    # Microbial
    println("saving figures")
    Plots.plot(tt1,At1,label="Av concentration profile when biculture",xaxis="Time(hr)",yaxis="Biomass(g/L)",framestyle=:box#=, legend=:topleft=#)
    Plots.plot!(tt1,St1,label="Se concentration profile when biculture")
    println(out_dir)
    Plots.savefig(out_dir * "Microbial_" * filename * "_X_Av0_$(@sprintf("%.2f",A0))_X_Se0_$(@sprintf("%.2f",S0))_N1_$(@sprintf("%.2f",s1)) N3_$(@sprintf("%.2f",s2)) N5_$(@sprintf("%.2f",s3)) N7_$(@sprintf("%.2f",s4)) N9_$(@sprintf("%.2f",s5)) N2_$(@sprintf("%.2f",a1)) N4_$(@sprintf("%.2f",a2)) N6_$(@sprintf("%.2f",a3)) N8_$(@sprintf("%.2f",a4)) N10_$(@sprintf("%.2f",a5)).pdf")
    Plots.plot(tt1,Ct1,label="Sucrose concentration profile when biculture",xaxis="Time(hr)",yaxis="Nutrient(g/L)",framestyle=:box#=, legend=:topleft=#)
    Plots.plot!(tt1,Nt1,label="Ammonia concentration profile when biculture")
    Plots.savefig(out_dir * "Nutrient_" * filename * "_X_Av0_$(@sprintf("%.2f",A0))_X_Se0_$(@sprintf("%.2f",S0))_N1_$(@sprintf("%.2f",s1)) N3_$(@sprintf("%.2f",s2)) N5_$(@sprintf("%.2f",s3)) N7_$(@sprintf("%.2f",s4)) N9_$(@sprintf("%.2f",s5)) N2_$(@sprintf("%.2f",a1)) N4_$(@sprintf("%.2f",a2)) N6_$(@sprintf("%.2f",a3)) N8_$(@sprintf("%.2f",a4)) N10_$(@sprintf("%.2f",a5)).pdf")

    # Store data into excel files
    println("writing plots to files")
    # println(out_dir)
    top_excel_file = out_dir * "\\" * filename * "_X_Av_$(@sprintf("%.2f",A0))_X_Se_$(@sprintf("%.2f",S0))_muS_$(@sprintf("%.2f",s1)) ksS_$(@sprintf("%.2f",s2)) ysxS_$(@sprintf("%.2f",s3)) yspS_$(@sprintf("%.2f",s4)) msS_$(@sprintf("%.2f",s5)) muA_$(@sprintf("%.2f",a1)) ksA_$(@sprintf("%.2f",a2)) ysxA_$(@sprintf("%.2f",a3)) yspA_$(@sprintf("%.2f",a4)) msA_$(@sprintf("%.2f",a5)).xlsx"
    column_names = ["times (hr)","Av","Se", "Sucrose", "Ammonia"]
    data=[tt1,At1,St1,Ct1,Nt1]
    # write to excel file
    XLSX.writetable(top_excel_file, data, column_names)
  
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

function Startup(A,S,N,C,tspan1,pS,pA) # batch
    # Monod Equation
    f(y,p,t)=[pA[1]*y[3]/(y[3]+pA[3])*y[1]-pA[2]*y[1], # X(Av)
            pS[1]*y[4]/(y[4]+pS[3])*y[2]-pS[2]*y[2], # X(Av)
            pS[1]*pS[5]*y[4]/(y[4]+pS[3])*y[2]-pA[1]*pA[4]*y[3]/(y[3]+pA[3])*y[1], # Sucrose
            pA[1]*pA[5]*y[3]/(y[3]+pA[3])*y[1]-pS[1]*pS[4]*y[4]/(y[4]+pS[3])*y[2]] # Ammonia
    
 

    prob=ODEProblem(f,[A,S,C,N],(0.0,tspan1))
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


