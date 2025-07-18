# Simulations of Dimensionless ODEs
# using DifferentialEquations, NLsolve, XLSX, Printf
using Plots, DifferentialEquations, NLsolve, XLSX, Printf

# function loadProcessData()
# function loadProcessData(N9,N10)
# function loadProcessData(TC)
function loadProcessData(N9,N10,N11,N12)
    global mu_maxA = 0.11
    global kdA = 0.2*0.11 # 5% of mu_maxA # N4
    global Ks_Av = 5*342.3*10^(-3) # 5mM of sucrose
    global Ys_Av = 1
    global Yp_Av = 2 # alpha_N
    # global mu_maxS = 0.0217
    # Assume mu_maxS = mu_maxA for the simplest case
    global mu_maxS = 0.11
    global mu_maxE = 0.2
    global kdS = kdA
    global kdE = 0.2*mu_maxE
    global Ks_Se = 5*17.03*10^(-3) # 5mM of ammonia
    global Ks_S_Ec = 1 
    global Ks_N_Ec = 0.3 
    global Ys_Se = 1  
    global Ys_S_Ec = 1  
    global Ys_N_Ec = 1  
    global Yp_Se = 0.5 # alpha_C
    global Yp_Ec = 0.3 # alpha_C
    global beta_Se = 0.2 
    global beta_Av = 0.2
    global beta_Ec = 0.3
    global D_Av = 0.00
    global D_Se = 0.00
    global D_Ec = 0.02
    global ms_Av = 0 # maintainace coefficient, we assume it is 0
    global mn_Se = 0 # maintainace coefficient, we assume it is 0
    global ms_Ec = 0 # maintainace coefficient, we assume it is 0
    global mn_Ec = 0 # maintainace coefficient, we assume it is 0

    V_se = 1
    V_av = 1
    V_ec = 1
    V = V_se+V_av+V_ec
    
    # Test only for coculture and follow the dimensionless number value assumption in DoE catch up meeting 1204_2 slides
    tc = 1/mu_maxS # N1 = 1
    # tc = TC # N1 = 1
    # C_S_c = Ks_Av
    # global C_S_c = Ks_Se
    global C_S_c = Ks_Av
    global C_NH4_c = Ks_Se
    global X_Se_c = Ks_Se*Ys_Se
    global X_Av_c = Ks_Av*Ys_Av
    global C_P_c = Yp_Ec*Ks_Av*Ys_S_Ec
    global X_Ec_c = Ks_Av*Ys_S_Ec
    N1 = tc*mu_maxS
    N2 = tc*mu_maxA
    N3 = tc*kdS
    N4 = tc*kdA
    N5 = Ks_Se/C_NH4_c
    N6 = Ks_Av/C_S_c
    # N7 = X_Se_c/C_NH4_c/Ys_Se
    # N8 = X_Av_c/C_S_c/Ys_Av
    # N9 = Yp_Se*X_Se_c/C_S_c
    # N10 = Yp_Av*X_Av_c/C_NH4_c
    N7 = 1
    N8 = 1
    # N9 = 1
    # N10 = 1
    # N11 = 0.1
    # N12 = 0.1
    # N11 = tc*beta_Se*X_Se_c/C_S_c
    # N12 = tc*beta_Av*X_Av_c/C_NH4_c
    # N11 = 0.00
    # N12 = 0.00

    N13 = tc*mu_maxE
    N14 = tc*kdE
    N15 = tc*D_Se
    N16 = tc*D_Av
    N17 = tc*D_Ec
    N18 = Ks_Se/C_NH4_c
    N19 = Ks_S_Ec/C_S_c
    N20 = 1
    N21 = Ks_Av*Ys_S_Ec/Ks_Se/Ys_N_Ec
    N22 = 1
    N23 = beta_Ec*tc/Yp_Ec
    N24 = V_se/V
    N25 = V_av/V
    N26 = V_ec/V
    N27 = tc*ms_Av*X_Av_c/C_S_c
    N28 = tc*ms_Ec*X_Ec_c/C_S_c
    N29 = tc*mn_Se*X_Se_c/C_NH4_c
    N30 = tc*mn_Ec*X_Ec_c/C_NH4_c

    # global A0 = 0.01
    # global S0 = 0.01
    global tspan1=20000
    global tspan2=100
    global out_dir="G:\\My Drive\\Research\\DOE project\\Modeling\\DimensionlessAnalysis\\test0422formanuscript\\"

    println("Parameters Loaded!")
    return [tc X_Se_c X_Av_c C_S_c C_NH4_c C_P_c X_Ec_c D_Ec D_Av], [N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12 N13 N14 N15 N16 N17 N18 N19 N20 N21 N22 N23 N24 N25 N26 N27 N28 N29 N30]
end

# function BicultureGrowth(N0, C0; filename = "NotSpecific") # Batch flow
# function BicultureGrowth(N0, C0, N9, N10; filename = "NotSpecific") # Batch flow
function BicultureGrowth(N0, C0, A0_given, S0_given, N9, N10, N11, N12; filename = "NotSpecific") # Batch flow
# function BicultureGrowth(N0, C0, A0_given, S0_given, TC; filename = "NotSpecific") # Batch flow
# function BicultureGrowth(N0, C0, A0_given, S0_given; filename = "BatchBiculture") # Batch flow
    # Critical_Numbers, Dimensionless_Numbers = loadProcessData()
    # Critical_Numbers, Dimensionless_Numbers = loadProcessData(N9,N10)
    # Critical_Numbers, Dimensionless_Numbers = loadProcessData(TC)
    Critical_Numbers, Dimensionless_Numbers = loadProcessData(N9,N10,N11,N12)
    println("Critical Numbers are: ", Critical_Numbers)
    
    Boundaries = Dimensionless_Numbers[9]*Dimensionless_Numbers[10]/Dimensionless_Numbers[7]/Dimensionless_Numbers[8]
    println("Boundaries = ", Boundaries)
    # if Boundaries != 1.0
    #     println("Round Dimensionless Numbers!")
    #     Dimensionless_Numbers[7] = round(X_Se_c/C_NH4_c/Ys_Se, digits = 2)
    #     Dimensionless_Numbers[8] = round(X_Av_c/C_S_c/Ys_Av, digits = 2)
    # end
    println("Dimensionless Numbers are: ", Dimensionless_Numbers)
    # global tt1,At1,St1,Ct1,Nt1=Startup1(A0,S0,N0,C0,tspan1,Critical_Numbers)
    # global tt1,At1,St1,Ct1,Nt1=Startup2(A0_given,S0_given,N0,C0,tspan1,Dimensionless_Numbers)
    global tt1,At1,St1,Ct1,Nt1,At2,St2=Startup3(A0_given,S0_given,N0,C0,tspan1,Dimensionless_Numbers) # At2 and St2 is the measured cell concentration
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
    c1 = Critical_Numbers[1]
    c2 = Critical_Numbers[2]
    c3 = Critical_Numbers[3]
    c4 = Critical_Numbers[4]
    c5 = Critical_Numbers[5]

    # Microbial
    println("saving figures")
    Plots.plot(tt1,At1,label="A",xaxis="Dimensionless Time",yaxis="Dimensionless Biomass Concentration",framestyle=:box#=, legend=:topleft=#)
    Plots.plot!(tt1,St1,label="B")
    println(out_dir)
    Plots.savefig(out_dir * "Microbial_" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.6f",C0))_N0_$(@sprintf("%.6f",N0))_N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).pdf")
    Plots.plot(tt1,Ct1,label="B's product",xaxis="Dimensionless Time",yaxis="Dimensionless Nutrient Concentration",framestyle=:box#=, legend=:topleft=#)
    Plots.plot!(tt1,Nt1,label="A's product")
    Plots.savefig(out_dir * "Nutrient_" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.6f",C0))_N0_$(@sprintf("%.6f",N0))_N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).pdf")

    plot(tt1,At1,label="A(B)",xaxis="Dimensionless Time",yaxis="Dimensionless Living Concentration", color=:blue)#, legend=:topright)
    plot!(twinx(),tt1, Nt1, yaxis="Dimensionless Nutrient Concentration", label="AP(BP)",color=:red)#), legend=:topleft)
    Plots.savefig(out_dir * "TEST_" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.6f",C0))_N0_$(@sprintf("%.6f",N0))_N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).pdf")


    # Store data into excel files
    println("writing plots to files")
    # println(out_dir)
    # top_excel_file = out_dir * "\\" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.4f",C0))_N0_$(@sprintf("%.4f",N0))_tc_$(@sprintf("%.2f",c1)) XSec_$(@sprintf("%.2f",c2)) XAvc_$(@sprintf("%.2f",c3)) CSc_$(@sprintf("%.2f",c4)) CAs_$(@sprintf("%.2f",c5)) YsA_$(@sprintf("%.2f",Ys_Av)) YsS_$(@sprintf("%.2f",Ys_Se)) YpA_$(@sprintf("%.2f",Yp_Av)) YpS_$(@sprintf("%.2f",Yp_Se)) N9_$(@sprintf("%.2f",Dimensionless_Numbers[9])) N10_$(@sprintf("%.2f",Dimensionless_Numbers[10])) N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).xlsx"
    top_excel_file = out_dir * "\\" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.6f",C0))_N0_$(@sprintf("%.6f",N0))_N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).xlsx"

    column_names = ["times (hr)","A(Av)","B(Se)", "BP(Sucrose)", "AP(Ammonia)", "A(Av measured)", "B(Se measured)"]
    data=[tt1,At1,St1,Ct1,Nt1, At2, St2]
    # write to excel file
    XLSX.writetable(top_excel_file, data, column_names)
    
    return tt1,At1,St1,Ct1,Nt1,At2, St2
  
end

function AllGrowth(N0,C0,P0,A0_given,S0_given,E0; filename = "Triculture") # Continuous flow
    tt1, At1, St1, Ct1, Nt1 = BicultureGrowth(N0, C0, A0_given, S0_given; filename = "NotSpecific")
    Critical_Numbers, Dimensionless_Numbers = loadProcessData()
    # Critical_Numbers, Dimensionless_Numbers = loadProcessData(N9,N10,N11,N12)
    c1 = Critical_Numbers[1]
    c2 = Critical_Numbers[2]
    c3 = Critical_Numbers[3]
    c4 = Critical_Numbers[4]
    c5 = Critical_Numbers[5]
    c8 = Critical_Numbers[8] # Dilution rate
    c9 = Critical_Numbers[9]
    # Incubate E.coli to start triculturing
    global tt,At, St, Et, Ct, Nt, Pt = Tripartite(E0,At1[end],St1[end],Ct1[end],Nt1[end],P0,tspan2,Dimensionless_Numbers)
    # E,A,S mean E.coli, Av, and Se
    # global dEdt=zeros(size(tt)[1])
    # global dAdt=zeros(size(tt)[1])
    # global dSdt=zeros(size(tt)[1])
    # global dCdt1=zeros(size(tt)[1])
    # global dCdt2=zeros(size(tt)[1])
    # global dCdt3=zeros(size(tt)[1])
    # global dNdt1=zeros(size(tt)[1])
    # global dNdt2=zeros(size(tt)[1])
    # global dNdt3=zeros(size(tt)[1])
    # global dPdt=zeros(size(tt)[1])
    size_t=size(tt)[1]+size(tt1)[1]
    # global final_dEdt=zeros(size_t)
    # global final_dAdt=zeros(size_t)
    # global final_dSdt=zeros(size_t)
    global final_Et=zeros(size_t)
    global final_At=zeros(size_t)
    global final_St=zeros(size_t)
    global final_Ct=zeros(size_t)
    global final_Nt=zeros(size_t)
    global final_Pt=zeros(size_t)

    
    ttt=cat(tt1,tt.+tspan1;dims=(1,1))
    # println(size(ttt)[1])
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

  
    final_Et=cat(zeros(size(tt1)[1]),Et;dims=(1,1))
    # plot(ttt,final_Et,title="E.coli concentration profile",xaxis="Time(hr)",yaxis="OD600",label=false)
    # savefig("Ecoli_test triculture.png")

    final_At=cat(At1,At;dims=(1,1))
    # plot(ttt,final_At,title="Av concentration profile",xaxis="Time(hr)",yaxis="Av(g/L)",label=false)
    # savefig("Av_test triculture.png")

    final_St=cat(St1,St;dims=(1,1))
    # plot(ttt,final_St,title="Se concentration profile",xaxis="Time(hr)",yaxis="Se(g/L)",label=false)
    # savefig("Se_test triculture.png")

    final_Ct=cat(Ct1,Ct;dims=(1,1))
    # plot(ttt,final_Ct,title="Sucrose concentration profile",xaxis="Time(hr)",yaxis="Sucrose(g/L)",label=false)
    # savefig("Sucrose_test triculture.png")

    final_Nt=cat(Nt1,Nt;dims=(1,1))
    # plot(ttt,final_Nt,title="Ammonia concentration profile",xaxis="Time(hr)",yaxis="Ammonia(g/L)",label=false)
    # savefig("Ammonia_test triculture.png")

    final_Pt=cat(zeros(size(tt1)[1]),Pt;dims=(1,1))
    # plot(ttt,final_Pt,title="Product concentration profile",xaxis="Time(hr)",yaxis="Isobutanol(g/L)",label=false)
    # savefig("Isobutanol_test.png")

    # Store data into excel files
    # println("writing plots to files")
    # top_excel_file = out_dir * "\\Profiles of All Microbial without inhibition.xlsx"
    # column_names = ["times (hr)","Ec","Av","Se", "Sucrose", "Ammonia","Isobutanol"]
    # data=[ttt,final_Et,final_At,final_St,final_Ct,final_Nt,final_Pt]
    # # write to excel file
    # XLSX.writetable(top_excel_file, data, column_names)

    # Microbial
    println("saving figures")
    Plots.plot(ttt,final_At,label="Av",xaxis="Dimensionless Time",yaxis="Dimensionless Biomass Concentration",framestyle=:box#=, legend=:topleft=#)
    Plots.plot!(ttt,final_St,label="Se")
    Plots.plot!(ttt,final_Et,label="Ec")
    println(out_dir)
    Plots.savefig(out_dir * "Final_Microbial_" * filename * "_t1_$(@sprintf("%.3f",tspan1))_S0_$(@sprintf("%.3f",C0))_tc_$(@sprintf("%.2f",c1)) XSec_$(@sprintf("%.2f",c2)) XAvc_$(@sprintf("%.2f",c3)) CSc_$(@sprintf("%.2f",c4)) CAs_$(@sprintf("%.2f",c5)) D_$(@sprintf("%.2f",c8))_$(@sprintf("%.2f",c9)) YsA_$(@sprintf("%.2f",Ys_Av)) YsS_$(@sprintf("%.2f",Ys_Se)) YpA_$(@sprintf("%.2f",Yp_Av)) YpS_$(@sprintf("%.2f",Yp_Se)) N9_$(@sprintf("%.2f",Dimensionless_Numbers[9])) N10_$(@sprintf("%.2f",Dimensionless_Numbers[10])) N11_$(@sprintf("%.3f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.3f",Dimensionless_Numbers[12])).pdf")
    Plots.plot(ttt,final_Ct,label="Sucrose",xaxis="Dimensionless Time",yaxis="Dimensionless Nutrient Concentration",framestyle=:box#=, legend=:topleft=#)
    Plots.plot!(ttt,final_Nt,label="Ammonia")
    Plots.plot!(ttt,final_Pt,label="Product")
    Plots.savefig(out_dir * "Final_Nutrient_" * filename * "_t1_$(@sprintf("%.3f",tspan1))_S0_$(@sprintf("%.3f",C0))_tc_$(@sprintf("%.2f",c1)) XSec_$(@sprintf("%.2f",c2)) XAvc_$(@sprintf("%.2f",c3)) CSc_$(@sprintf("%.2f",c4)) CAs_$(@sprintf("%.2f",c5)) D_$(@sprintf("%.2f",c8))_$(@sprintf("%.2f",c9)) YsA_$(@sprintf("%.2f",Ys_Av)) YsS_$(@sprintf("%.2f",Ys_Se)) YpA_$(@sprintf("%.2f",Yp_Av)) YpS_$(@sprintf("%.2f",Yp_Se)) N9_$(@sprintf("%.2f",Dimensionless_Numbers[9])) N10_$(@sprintf("%.2f",Dimensionless_Numbers[10])) N11_$(@sprintf("%.3f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.3f",Dimensionless_Numbers[12])).pdf")

    # Store data into excel files
    println("writing plots to files")
    # println(out_dir)
    top_excel_file = out_dir * "\\" * filename * "_t1_$(@sprintf("%.3f",tspan1))_S0_$(@sprintf("%.3f",C0))_tc_$(@sprintf("%.2f",c1)) XSec_$(@sprintf("%.2f",c2)) XAvc_$(@sprintf("%.2f",c3)) CSc_$(@sprintf("%.2f",c4)) CAs_$(@sprintf("%.2f",c5)) D_$(@sprintf("%.2f",c8))_$(@sprintf("%.2f",c9)) YsA_$(@sprintf("%.2f",Ys_Av)) YsS_$(@sprintf("%.2f",Ys_Se)) YpA_$(@sprintf("%.2f",Yp_Av)) YpS_$(@sprintf("%.2f",Yp_Se)) N9_$(@sprintf("%.2f",Dimensionless_Numbers[9])) N10_$(@sprintf("%.2f",Dimensionless_Numbers[10])) N11_$(@sprintf("%.3f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.3f",Dimensionless_Numbers[12])).xlsx"
    column_names = ["times (hr)","Av","Se", "Ec", "Sucrose", "Ammonia", "Product"]
    data=[ttt,final_At, final_St, final_Et, final_Ct, final_Nt, final_Pt]
    # write to excel file
    XLSX.writetable(top_excel_file, data, column_names)
  
end

function Startup1(A,S,N,C,tspan1,cr) # batch
    # Monod Equation
    # [tc X_Se_c X_Av_c C_S_c C_NH4_c]
    # Don't use it because it has numerical error: e.g. when X_Se_c=C_S_c, X_Av_c=C_NH4_c, N7*N8=0.999999, which leads to N9*N10/N7/N8>1
    f(y,p,t)=[cr[1]*mu_maxA*y[3]/(y[3]+Ks_Av/cr[4])*y[1] - cr[1]*kdA*y[1], # X(Av)
              cr[1]*mu_maxS*y[4]/(y[4]+Ks_Se/cr[5])*y[2] - cr[1]*kdS*y[2], # X(Se)
              cr[1]*mu_maxS*Yp_Se*cr[2]/cr[4]*y[4]/(y[4]+Ks_Se/cr[5])*y[2] + beta_Se*cr[2]*cr[1]/cr[4]*y[2] - cr[1]*mu_maxA*cr[3]/cr[4]/Ys_Av*y[3]/(y[3]+Ks_Av/cr[4])*y[1], # Sucrose
              cr[1]*mu_maxA*Yp_Av*cr[3]/cr[5]*y[3]/(y[3]+Ks_Av/cr[4])*y[1] + beta_Av*cr[3]*cr[1]/cr[5]*y[1] - cr[1]*mu_maxS*cr[2]/cr[5]/Ys_Se*y[4]/(y[4]+Ks_Se/cr[5])*y[2]] # Ammonia    
 

    prob=ODEProblem(f,[A,S,C,N],(0.0,tspan1))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    return a,A[1,:],A[2,:],A[3,:],A[4,:]
    # At,St,Ct,Nt
end

function Startup2(A,S,N,C,tspan1,dn) # batch
    # Monod Equation
    # [tc X_Se_c X_Av_c C_S_c C_NH4_c]
    f(y,p,t)=[dn[2]*y[3]/(y[3]+dn[6])*y[1] - dn[4]*y[1], # X(Av)
              dn[1]*y[4]/(y[4]+dn[5])*y[2] - dn[3]*y[2], # X(Se)
              dn[1]*dn[9]*y[4]/(y[4]+dn[5])*y[2] - dn[2]*dn[8]*y[3]/(y[3]+dn[6])*y[1], # Sucrose
              dn[2]*dn[10]*y[3]/(y[3]+dn[6])*y[1] - dn[1]*dn[7]*y[4]/(y[4]+dn[5])*y[2]] # Ammonia
    

    prob=ODEProblem(f,[A,S,C,N],(0.0,tspan1))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    return a,A[1,:],A[2,:],A[3,:],A[4,:]
    # At,St,Ct,Nt
end

function Startup3(A,S,N,C,tspan1,dn) # batch
    # Monod Equation
    # Add nongrowth related term
    # [tc X_Se_c X_Av_c C_S_c C_NH4_c]
    f(y,p,t)=[dn[2]*y[3]/(y[3]+dn[6])*y[1] - dn[4]*y[1], # X(Av)
              dn[1]*y[4]/(y[4]+dn[5])*y[2] - dn[3]*y[2], # X(Se)
              dn[11]*y[2]- dn[2]*dn[8]*y[3]/(y[3]+dn[6])*y[1], # Sucrose
              dn[12]*y[1] - dn[1]*dn[7]*y[4]/(y[4]+dn[5])*y[2], # Ammonia
              dn[2]*y[3]/(y[3]+dn[6])*y[1], # measured X(Av) the total amount of dry biomass including live and dead cells
              dn[1]*y[4]/(y[4]+dn[5])*y[2]] # measured X(Se) 
    

    prob=ODEProblem(f,[A,S,C,N,A,S],(0.0,tspan1))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    return a,A[1,:],A[2,:],A[3,:],A[4,:], A[5,:],A[6,:]
    # At,St,Ct,Nt
end

function Tripartite(E,A,S,C,N,P,tspan2,dn) # Use one ODE solver to solve the whole system
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
    # f(y,p,t)=[(mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1]-D[1]*y[1], # X(E.coli)
    #      (mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]-D[2]*y[2],# X(Av)
    #      (mu_maxS*y[5]/(ksS+y[5]) - kdS)*y[3]-D[3]*y[3],# X(Se)
    #      1/Vt*max(y[4],0)/y[4]*(-(mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE + msE)*y[1] - ((mu_maxA*y[4]/(ksA+y[4]) - kdA)/ysxA+msA)*y[2] + yspS*(mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS*y[3]  - sum(D)*y[4]), # Sucrose
    #      1/Vt*max(y[5],0)/y[5]*(yspA*(mu_maxA*y[4]/(ksA+y[4]) - kdA)*y[2]/ysxA - (mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])/ysxE_NH4 + msE_NH4)*y[1] - ((mu_maxS*y[5]/(ksS+y[5]) - kdS)/ysxS+msS)*y[3] - sum(D)*y[5]), # Ammonia
    #      1/Vt*(max((mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/((mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5]) - kdE)*y[1])*mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])*ysp_g/ysxE + (1-max((mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1],0)/((mu_maxE*y[4]/(ksE+y[4])*y[5]/(ksE_NH4+y[5])-kdE)*y[1]))*ysp_m*msE)*y[1] - sum(D)*y[6]] # Product

    # Dimensionless ODEs
    # two substrate limited, without inhibition, no maintainace coefficient
    f(y,p,t)=[(dn[2]*y[4]/(y[4]+dn[6])-dn[4]-dn[16])*y[1], # X(Av)
              (dn[1]*y[5]/(y[5]+dn[5])-dn[3] - dn[15])*y[2], # X(Se)
              (dn[13]*y[4]/(y[4]+dn[19])*y[5]/(y[5]+dn[18])-dn[14]-dn[17])*y[3], # X(Ec)
              dn[24]*max(dn[1]*dn[9]*y[5]/(y[5]+dn[5])-dn[3]*dn[9]-dn[15]*dn[9]+dn[11],0)*y[2] - dn[25]*dn[2]*dn[8]*max(y[4],0)/(y[4]+dn[6])*y[1] - dn[27]*dn[25]*y[1]*max(y[4],0)/y[4] - dn[26]*dn[13]*dn[20]*max(y[4],0)/(y[4]+dn[19])*max(y[5],0)/(y[5]+dn[18])*y[3] - dn[26]*dn[28]*y[3]*max(y[4],0)/y[4] - max(y[4],0)*(dn[24]*dn[15]+dn[25]*dn[16]+dn[26]*dn[17]), # Sucrose
              dn[25]*max(dn[10]*dn[2]*y[4]/(y[4]+dn[6])-dn[10]*dn[2]-dn[16]*dn[10]+dn[12],0)*y[1] - dn[24]*dn[1]*dn[7]*max(y[5],0)/(y[5]+dn[5])*y[2] - dn[29]*dn[24]*y[2]*max(y[5],0)/y[5] - dn[26]*dn[13]*dn[21]*max(y[4],0)/(y[4]+dn[19])*max(y[5],0)/(y[5]+dn[18])*y[3] - dn[26]*dn[30]*y[3]*max(y[5],0)/y[5] - max(y[5],0)*(dn[24]*dn[15]+dn[25]*dn[16]+dn[26]*dn[17]), # Ammonia
              dn[26]*max(dn[22]*dn[13]*y[4]/(y[4]+dn[19])*y[5]/(y[5]+dn[18])-dn[22]*dn[14]-dn[22]*dn[17]+dn[23],0)*max(y[3],0) - max(y[6],0)*(dn[24]*dn[15]+dn[25]*dn[16]+dn[26]*dn[17]) ]# Product 

    prob=ODEProblem(f,[A,S,E,C,N,P],(0.0,tspan2))
    # PositiveDomain(S=nothing;save=true,abstol=nothing,scalefactor=nothing)
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    a=soln.t
    A=Array(soln)
    # plot(a,A[1,:],title="Microbial concentration profile",xaxis="Time(hr)",yaxis="X(g/L)",label=false)
    # plot(a,A[2,:],title="Substrate concentration profile",xaxis="Time(hr)",yaxis="S(g/L)",label=false)
    # plot(a,A[3,:],title="Product concentration profile",xaxis="Time(hr)",yaxis="P(g/L)",label=false)
    return a,A[1,:],A[2,:],A[3,:],A[4,:],A[5,:],A[6,:]
    # Et,At,St,Pt,Ct,Nt
end



# BicultureGrowth(0.2, 0.2, 0.01, 0.01, 0, 0, 0.21, 0.21; filename = "Nongrowthonly_P_steady_0.25_")

# tt1 = [0, 8.944799847, 98.39279831, 992.872783, 2000]
# At1 = [0.01, 0.01, 0.01, 0.01, 0.01]
# St1 = [0.01, 0.01, 0.01, 0.01, 0.01]
# Ct1 = [0.25, 0.25, 0.25, 0.25, 0.25]
# Nt1 = [0.25, 0.25, 0.25, 0.25, 0.25]

# Plots.plot(tt1,At1,label="A",xaxis="Dimensionless Time",yaxis="Dimensionless Biomass Concentration",framestyle=:box#=, legend=:topleft=#, ylims=[0.0,0.02])
# Plots.plot!(tt1,St1,label="B")

# # Plots.savefig(out_dir * "Microbial_" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.4f",C0))_N0_$(@sprintf("%.4f",N0))_tc_$(@sprintf("%.2f",c1)) XSec_$(@sprintf("%.2f",c2)) XAvc_$(@sprintf("%.2f",c3)) CSc_$(@sprintf("%.2f",c4)) CAs_$(@sprintf("%.2f",c5)) YsA_$(@sprintf("%.2f",Ys_Av)) YsS_$(@sprintf("%.2f",Ys_Se)) YpA_$(@sprintf("%.2f",Yp_Av)) YpS_$(@sprintf("%.2f",Yp_Se)) N9_$(@sprintf("%.2f",Dimensionless_Numbers[9])) N10_$(@sprintf("%.2f",Dimensionless_Numbers[10])) N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).pdf")
# Plots.savefig(out_dir * "Microbial_" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.6f",C0))_N0_$(@sprintf("%.6f",N0))_N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).pdf")
# Plots.plot(tt1,Ct1,label="B's product",xaxis="Dimensionless Time",yaxis="Dimensionless Nutrient Concentration",framestyle=:box#=, legend=:topleft=#, ylims=[0.15,0.35])
# Plots.plot!(tt1,Nt1,label="A's product")
# # Plots.savefig(out_dir * "Nutrient_" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.4f",C0))_N0_$(@sprintf("%.4f",N0))_tc_$(@sprintf("%.2f",c1)) XSec_$(@sprintf("%.2f",c2)) XAvc_$(@sprintf("%.2f",c3)) CSc_$(@sprintf("%.2f",c4)) CAs_$(@sprintf("%.2f",c5)) YsA_$(@sprintf("%.2f",Ys_Av)) YsS_$(@sprintf("%.2f",Ys_Se)) YpA_$(@sprintf("%.2f",Yp_Av)) YpS_$(@sprintf("%.2f",Yp_Se)) N9_$(@sprintf("%.2f",Dimensionless_Numbers[9])) N10_$(@sprintf("%.2f",Dimensionless_Numbers[10])) N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).pdf")
# Plots.savefig(out_dir * "Nutrient_" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.6f",C0))_N0_$(@sprintf("%.6f",N0))_N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).pdf")

# 03/31/2025
# BicultureGrowth(0.25, 0.25, 0.01, 0.01, 0, 0, 0.1, 0.4; filename = "Nongrowthonly_N11_0.1_N12_0.4_")


# CA_0 = LinRange(0.1, 1, 10) # Ammonia
# CB_0 = LinRange(0.1, 1, 10) # Sucrose
# for k in eachindex(CA_0)
#     for l in eachindex(CB_0)
#         BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.2, 0.2; filename = "Nongrowthonly_N11N12equal0.04_")
#     end
# end

# BicultureGrowth(1/9, 0.25, 0.01, 0.01, 0, 0, 0.1, 0.1; filename = "Nongrowthonly_N11N12_0.01_")

# CA_0 = LinRange(0.1, 1, 10) # Ammonia
# CB_0 = LinRange(0.1, 1, 10) # Sucrose
# for k in eachindex(CA_0)
#     for l in eachindex(CB_0)
#         BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.1, 0.1; filename = "Nongrowthonly_N11N12_0.01_")
#     end
# end

# # N11=N12=0.2
# sum_PA_0_PB_0_1 = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
# xA_steady_1 = [0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46]

# # N11=N12=0.1
# sum_PA_0_PB_0_2 = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
# xA_steady_2 = [0.0294, 0.0794, 0.1294, 0.1794, 0.2294, 0.2794, 0.3294, 0.3794, 0.4294, 0.4794, 0.5294]

# # N11=N12=0.4
# sum_PA_0_PB_0_3 = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
# xA_steady_3 = [0.051666601, 0.101666665, 0.151666632, 0.201666794, 0.251666858, 0.301666922, 0.351666986]

# plot(sum_PA_0_PB_0_1, xA_steady_1, xlabel="PA + PB intial conc. (dimensionless)", ylabel="XA(XB) cell conc. (dimensionless)", label="N11=0.2, N12=0.2, N3=N4=0.2")
# plot!(sum_PA_0_PB_0_2, xA_steady_2, label="N11=0.1, N12=0.1, N3=N4=0.1")
# plot!(sum_PA_0_PB_0_3, xA_steady_3, label="N11=0.4, N12=0.4, N3=N4=0.4")

# # 04/03/2025
# # BicultureGrowth(0.25, 0.25, 0.01, 0.01, 0, 0, 0.1, 0.4; filename = "Nongrowthonly_N11_0.1_N12_0.4_")
# CA_0 = LinRange(0.1, 1, 10) # Ammonia
# CB_0 = LinRange(0.1, 1, 10) # Sucrose
# for k in eachindex(CA_0)
#     for l in eachindex(CB_0)
#         BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.4, 0.1; filename = "Nongrowthonly_N11_0.4_N12_0.1_")
#     end
# end

# # 04/06/2025
# BicultureGrowth(0.25, 9/11, 1, 1, 0, 0, 0.3, 0.3; filename = "Nongrowthonly_N11_0.3_N12_0.3_N4_0.45")

# CA_0 = LinRange(0.1, 1, 10) # Ammonia
# CB_0 = LinRange(0.3, 3, 10) # Sucrose
# for k in eachindex(CA_0)
#     for l in eachindex(CB_0)
#         BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.2, 0.2; filename = "Nongrowthonly_N11_0.2_N12_0.2_N3_0.4_N4_0.1")
#     end
# end

# BicultureGrowth(0.25, 1/19, 0.01, 0.01, 0, 0, 0.2, 0.05; filename = "Nongrowthonly_N11_0.2_N12_0.05_N4_0.05")
# BicultureGrowth(0.25, 1/19, 0.01, 0.01, 0, 0, 0.05, 0.2; filename = "Nongrowthonly_N11_0.05_N12_0.2_N4_0.05")

# BicultureGrowth(0.25, 0.128, 0.01, 0.01, 0, 0, 0.15, 0.15; filename = "Nongrowthonly_N11_0.15_N12_0.15_N4_0.1125")

# # 04/09/2025
# # N11 = 0.25, N12 = 0.16
# CA_0 = LinRange(0.1, 1, 10) # Ammonia
# CB_0 = LinRange(0.1, 1, 10) # Sucrose
# for k in eachindex(CA_0)
#     for l in eachindex(CB_0)
#         BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.25, 0.16; filename = "Nongrowthonly_N11_0.25_N12_0.16_N3_0.2_N4_0.2")
#     end
# end

# # N11 = 0.3, N12 = 4/30
# CA_0 = LinRange(0.1, 1, 10) # Ammonia
# CB_0 = LinRange(0.1, 1, 10) # Sucrose
# for k in eachindex(CA_0)
#     for l in eachindex(CB_0)
#         BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.3, 0.133; filename = "Nongrowthonly_N11_0.3_N12_0.133_N3_0.2_N4_0.2")
#     end
# end

# # N3 = 0.16, N4 = 0.25
# CA_0 = LinRange(0.1, 1, 10) # Ammonia
# CB_0 = LinRange(0.1, 1, 10) # Sucrose
# for k in eachindex(CA_0)
#     for l in eachindex(CB_0)
#         BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.2, 0.2; filename = "Nongrowthonly_N11_0.2_N12_0.2_N3_0.16_N4_0.25")
#     end
# end

# # N3 = 4/30, N4 = 0.3
# CA_0 = LinRange(0.1, 1, 10) # Ammonia
# CB_0 = LinRange(0.2, 2, 10) # Sucrose
# for k in eachindex(CA_0)
#     for l in eachindex(CB_0)
#         BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.2, 0.2; filename = "Nongrowthonly_N11_0.2_N12_0.2_N3_0.133_N4_0.3")
#     end
# end

# 04/10/2025
# N3 = 0.4, N4 = 0.4, N11=0.4, N12=0.4
CA_0 = LinRange(0.5, 5, 10) # Ammonia
CB_0 = LinRange(0.5, 5, 10) # Sucrose
for k in eachindex(CA_0)
    for l in eachindex(CB_0)
        BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.4, 0.4; filename = "Nongrowthonly_N11_0.4_N12_0.4_N3_0.4_N4_0.4")
    end
end

# N3 = 0.4, N4 = 0.4, N11=0.2, N12=0.8
CA_0 = LinRange(0.5, 5, 10) # Ammonia
CB_0 = LinRange(0.5, 5, 10) # Sucrose
for k in eachindex(CA_0)
    for l in eachindex(CB_0)
        BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.2, 0.8; filename = "Nongrowthonly_N11_0.2_N12_0.8_N3_0.4_N4_0.4")
    end
end

# N3 = 0.4, N4 = 0.4, N11=0.1, N12=1.6
CA_0 = LinRange(0.5, 5, 10) # Ammonia
CB_0 = LinRange(0.5, 5, 10) # Sucrose
for k in eachindex(CA_0)
    for l in eachindex(CB_0)
        BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.1, 1.6; filename = "Nongrowthonly_N11_0.1_N12_1.6_N3_0.4_N4_0.4")
    end
end

# N3 = 0.2, N4 = 0.4, N11=0.2, N12=0.4
CA_0 = LinRange(0.1, 1, 10) # Ammonia
CB_0 = LinRange(0.5, 5, 10) # Sucrose
for k in eachindex(CA_0)
    for l in eachindex(CB_0)
        BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.2, 0.4; filename = "Nongrowthonly_N11_0.2_N12_0.4_N3_0.2_N4_0.4")
    end
end

# N3 = 0.2, N4 = 0.1, N11=0.2, N12=0.1
CA_0 = LinRange(0.1, 1, 10) # Ammonia
CB_0 = LinRange(0.1, 1, 10) # Sucrose
for k in eachindex(CA_0)
    for l in eachindex(CB_0)
        BicultureGrowth(CA_0[k], CB_0[l], 0.01, 0.01, 0, 0, 0.2, 0.1; filename = "Nongrowthonly_N11_0.2_N12_0.1_N3_0.2_N4_0.1")
    end
end

BicultureGrowth(0.25, 0.25, 0.01, 0.01, 0, 0, 0.25, 0.25; filename = "Nongrowthonly_N11_0.25_N12_0.25_N3_0.2_N4_0.2")
BicultureGrowth(0.25, 0.25, 0.01, 0.01, 0, 0, 0.15, 0.15; filename = "Nongrowthonly_N11_0.15_N12_0.15_N3_0.2_N4_0.2")

# 04/18/2025
BicultureGrowth(0.1, 0.1, 0.01, 0.01, 0, 0, 0.21, 0.21; filename = "Nongrowthonly_N11_0.21_N12_0.21_N3_0.2_N4_0.2")

# 04/23/2025
BicultureGrowth(0.2255, 0.2255, 0.01, 0.01, 0, 0, 0.21, 0.21; filename = "Nongrowthonly_20000h_N11_0.21_N12_0.21_N3_0.2_N4_0.2")
