using Plots, Printf, XLSX
# using DataFrames
using Plots.PlotMeasures

xf = XLSX.readxlsx("G:\\My Drive\\Research\\DOE project\\Modeling\\coculture\\Monod equation\\Parameter permutation\\Excelfiles\\muS_0.11 ksS_0.11 ysxS_0.20 yspS_0.20 msS_0.11 muA_0.11 ksA_0.02 ysxA_0.11 yspA_0.02 msA_0.02.xlsx")
function MyPlots(xf)
    sh = xf["Sheet1"]
    num_rows, num_cols = size(sh[:])
    # println(string(num_rows) * " " * string(num_cols))
    range = "G3:G" * string(num_rows)
    # println(range)
    t = zeros(num_rows)
    A = zeros(num_rows)
    S = zeros(num_rows)
    C = zeros(num_rows)
    N = zeros(num_rows)

    t = convert(Array{Float64,2},sh["A3:A" * string(num_rows)]) 
    A = convert(Array{Float64,2},sh["B3:B" * string(num_rows)]) 
    S = convert(Array{Float64,2},sh["C3:C" * string(num_rows)]) 
    C = convert(Array{Float64,2},sh["D3:D" * string(num_rows)]) 
    N = convert(Array{Float64,2},sh["E3:E" * string(num_rows)]) 
    
    p=plot(t,A,label="A.v.",xlabel="times (h)",ylabel="Concentration (g/L)",legend=:topleft,framestyle=:box)
    plot!(t,S, label="S.e.")
    # savefig("Av and Se.pdf")
    plot!(t,C/10, label="Sucrose/10")
    plot!(t,N/10, label="Ammonia/10")
    
    display(p)
    # savefig("All profile.pdf")
    # p2=plot(tt,E,label="E.coli",xlabel="times (h)",ylabel="Concentration (g/L)",legend=:topleft,framestyle=:box)
    # # plot!(tt,E,label="E.coli")
    # plot!(tt,P, label="Isobutanol")
    # savefig("With inhibition only E and P up to 40 hrs.pdf")
    # p2=plot(tt,A, label="A.v.",xlabel="times (h)",ylabel="Concentration (g)",legend=:topright,framestyle=:box)
    # # plot!(tt,A, label="A.v.")
    # plot!(tt,S, label="S.e.")
    # # p3=plot(tt,C, label="Sucrose",xlabel="times (h)",ylabel="Concentration (g)",legend=:topright,framestyle=:box)
    # plot!(tt,C, label="Sucrose")
    # plot!(tt,N, label="Ammonia")
        
    # p_all=plot(p1,p2,p3,p4,layout=(2,2),legend=:topright,xtickfontsize=6,ytickfontsize=6,xguidefontsize=8,yguidefontsize=8,framestyle=:box)
    # display(p2)
    # savefig("Without inhibition.pdf")
    # savefig("With inhibition.pdf")
end

MyPlots(xf)
