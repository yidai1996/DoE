using Plots, Printf, XLSX
# using DataFrames
# using Plots.PlotMeasures

# xf = XLSX.readxlsx("G:\\My Drive\\Research\\DOE project\\Modeling\\coculture\\Monod equation\\Parameter permutation\\Excelfiles\\muS_0.11 ksS_0.11 ysxS_0.20 yspS_0.20 msS_0.11 muA_0.11 ksA_0.02 ysxA_0.11 yspA_0.02 msA_0.02.xlsx")
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

# MyPlots(xf)

using PyPlot, XLSX
xf = XLSX.readxlsx("G:\\My Drive\\Research\\DOE project\\Modeling\\3D plots.xlsx")

function D3Plots(xf)
    sh = xf["Sheet1"]
    num_rows, num_cols = size(sh[:])
    println(num_rows)
    # println(string(num_rows) * " " * string(num_cols))
    # range = "G2:G" * string(num_rows)
    # println(range)
    C = zeros(num_rows-1)
    N = zeros(num_rows-1)
    N9N10 = zeros(num_rows-1)

    C = convert(Array{Float64,2},sh["A2:A" * string(num_rows-1)]) 
    N = convert(Array{Float64,2},sh["B2:B" * string(num_rows-1)]) 
    N9N10 = convert(Array{Float64,2},sh["E2:E" * string(num_rows-1)]) 

    surf(C,N,N9N10)
    # p = PyPlot.plot(C,N,N9N10,st=:surface,camera=(-30,30))
    # savefig("G:\\My Drive\\Research\\DOE project\\Modeling\\3D plots test.pdf")
    # display(p)

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
D3Plots(xf)

function Death_impact()
    x = range(0.05, 0.5, length=11)
    y=x./(1 .-x)
    plot(x,y,xlabel="death rate/max growth rate", ylabel="steady state of nutrient concentration", title="Death rate influence (N9N10=1)")
    savefig("G:\\My Drive\\Research\\DOE project\\Modeling\\DimensionlessAnalysis\\test0129deathrate\\Death_Cs_figure.pdf")
    savefig("G:\\My Drive\\Research\\DOE project\\Modeling\\DimensionlessAnalysis\\test0129deathrate\\Death_Cs_figure.png")
end
Death_impact()

# parameter space
using Plots
x = [1, 1, 1.1, 1.2, 1.3, 2.0]
y = [1,0.5, 0.45, 0.43, 0.41, 0.38]

# Manuscript paper
using XLSX
using Plots

# 三个独立的 Excel 文件路径
file_directory = "G:\\My Drive\\Research\\DOE project\\Modeling\\DimensionlessAnalysis\\manuscriptplots\\N9N10_1\\"
file_paths = ["G:\\My Drive\\Research\\DOE project\\Modeling\\DimensionlessAnalysis\\manuscriptplots\\N9N10_1\\Hybrid__XA0_0.010_XB0_0.010_S0_0.240000_N0_0.240000_N9_1.00_N10_1.00_N11_0.000 N12_0.000.xlsx", 
                "G:\\My Drive\\Research\\DOE project\\Modeling\\DimensionlessAnalysis\\manuscriptplots\\N9N10_1\\Hybrid__XA0_0.010_XB0_0.010_S0_0.250000_N0_0.250000_N9_1.00_N10_1.00_N11_0.000 N12_0.000.xlsx", 
                "G:\\My Drive\\Research\\DOE project\\Modeling\\DimensionlessAnalysis\\manuscriptplots\\N9N10_1\\Hybrid__XA0_0.010_XB0_0.010_S0_0.260000_N0_0.260000_N9_1.00_N10_1.00_N11_0.000 N12_0.000.xlsx"]

# 用来存储各子图
plot_list = Vector{Plots.Plot{Plots.GRBackend}}()

    # 打开工作簿并自动选取第一个工作表
    xf = XLSX.readxlsx(file_paths[1])
    sh = xf["Sheet1"]

    # 假设第一列是 x，第二列是 y1，第三列是 y2，第一行是表头
    x = convert(Array{Float64,2},sh["A2:A22"]) 
    y1 = convert(Array{Float64,2},sh["B2:B22"]) # BIOMASS 
    y2 = convert(Array{Float64,2},sh["D2:D22"]) # NUTRIENT 

    # 创建主图（左侧 Y 轴）
    # 1) 在左轴绘制 y1
    p = plot(
        x, y1;
        label = "Y1",
        xlabel = "X",
        ylabel = "Y1",
        lw = 2,
        axislinewidth = 2,
        framestyle = :box,
        legend = :topleft
    )
    # 2) 在同一图上，使用 yaxis=:right 绘制 y2，并设置其范围
    plot!(
        p,
        x, y2;
        label = "Y2",
        yaxis = :right,
        ylabel = "Y2",
        ylim = (0.2, 0.3),
        lw = 2,
        axislinewidth = 2,
        legend = :topright
    )

    push!(plot_list, p)


# 将三幅子图按 (3, 1) 布局合并
final_plot = plot(plot_list..., layout = (3, 1), size = (800, 1200))

display(final_plot)

plot(tt1,At1,label="A(B)",xaxis=L"N_{9} \times N_{10}",yaxis="Dimensionless Living Concentration", color=:blue)#, legend=:topright)
plot!(twinx(),tt1, Nt1, yaxis="Dimensionless Nutrient Concentration", label="AP(BP)",color=:red)#), legend=:topleft)
Plots.savefig(out_dir * "TEST_" * filename * "_XA0_$(@sprintf("%.3f",A0_given))_XB0_$(@sprintf("%.3f",S0_given))_S0_$(@sprintf("%.6f",C0))_N0_$(@sprintf("%.6f",N0))_N9_$(@sprintf("%.2f",Dimensionless_Numbers[9]))_N10_$(@sprintf("%.2f",Dimensionless_Numbers[10]))_N11_$(@sprintf("%.2f",Dimensionless_Numbers[11])) N12_$(@sprintf("%.2f",Dimensionless_Numbers[12])).pdf")

using LaTeXStrings
plot([1.0,1.0], [0.25, 0.5], xaxis = L"N_{9} \times N_{10}", yaxis="Initial Conc. of Nutrient (Dimensionless)", color=:coral, ylims=[0, 0.5], xlims=[0.5,2.0], label="",  lw = 2)
plot!([1.0, 1.15, 1.3, 1.5, 2.0], [0.25, 0.22, 0.21, 0.2, 0.195], label="", color=:blue, lw=2)
savefig("parameterregion_growthonly.png")

plot([1.0,1.0], [0.25, 0.5], xaxis = L"N_{11} \times N_{12}", yaxis="Initial Conc. of Nutrient (Dimensionless)", color=:coral, ylims=[0, 0.5], xlims=[0.5,2.0], label="",  lw = 2)
plot!([1.0, 1.15, 1.3, 1.5, 2.0], [0.25, 0.23, 0.22, 0.215, 0.21], label="", color=:coral, lw=2)
savefig("parameterregion_momgrowthonly.png")