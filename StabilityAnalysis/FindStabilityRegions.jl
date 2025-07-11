using XLSX, DataFrames, CSV

# Distinguish stability by checking the last two points in all growth curves
out_dir = "C:\\Users\\yid\\TemporaryResearchDataStorage\\doe29_tspan_200\\excel\\"
index_files = readdir(out_dir);
# build a dataframe for C/N and stability (growth curve regions)
df1 = DataFrame([Float64[], Float64[], Float64[]], ["Sucrose","Ammonia","Stability"])
df2 = DataFrame([Float64[], Float64[], Float64[]], ["Sucrose","Ammonia","Stability"])
count = 0
for i in eachindex(index_files)
    count += 1
    # search all excel files
    # find the last two points of biomass concentration
    ex = XLSX.readxlsx(out_dir * index_files[i])
    sh = ex["Sheet1"]
    num_rows, num_cols = size(sh[:])
    t = sh["A2:A"*string(num_rows)]
    t_checkpoint_row = Tuple(findall(x -> x>=50, t)[1])[1]
    Av = sh["B2:B"*string(num_rows)]
    N = sh["E2:E"*string(num_rows)]
    C = sh["D2:D"*string(num_rows)]
    Av1 = sh["B" * string(num_rows-1)]
    Av2 = sh["B" * string(num_rows)]
    Se1 = sh["C" * string(num_rows-1)]
    Se2 = sh["C" * string(num_rows)]
    N_start_1 = sh["E2"]
    N_start_2 = sh["E3"]
    N_end = N[end]
    C_end = C[end]

    # println("Av=",[Av1 Av2])
    # println("Se=",[Se1 Se2])
    # Check the stability
    check_exponential = (Av2-Av1)/(t[end]-t[end-1])
    println(check_exponential)
    if check_exponential > 1
        # going to exponential growth
        push!(df1,[sh["D2"], sh["E2"], 0])
    elseif (Av[t_checkpoint_row]-Av[t_checkpoint_row-1])/(t[t_checkpoint_row]-t[t_checkpoint_row-1])<0
        # has a bump for Av growth
        push!(df1,[sh["D2"], sh["E2"], 1])
    else # no bumps for Av growth and going to stationary stage
        push!(df1,[sh["D2"], sh["E2"], 2])
    end

    # For only Monod model
    if N_end-N[end-1]<0.001 
        push!(df2,[sh["D2"], sh["E2"], 0])
        
    else # not stable, C_N doesn't go to the stable point
        push!(df2,[sh["D2"], sh["E2"], 1])
    end

    # # For monod-LP model
    # if (N_end-N[end-1])/(t[end]-t[end-1]) > 1
    #     if C_end < 0.03
    #         # C is stable 
    #         push!(df2,[sh["D2"], sh["E2"], 0])
    #     else # c NOT TOUCH THE STABLE FINAL STATE
    #         println(C_end)
    #         push!(df2,[sh["D2"], sh["E2"], 1])
    #     end
    # else # N is not increasing exponentially
    #     push!(df2,[sh["D2"], sh["E2"], 2])
    # end
    println(i)
end
println(count)
println(df1)
println(df2)
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\doe29_tspan_200\\StabilityAnalysisResults\\StabilityRegionsMicrobial.csv",df1)
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\doe29_tspan_200\\StabilityAnalysisResults\\StabilityRegionsNutrient.csv",df2)


# Plot the stability regions
using Plots
# using StatsPlots
# read data from CSV file by making a dataframe
df_nutrient = CSV.read("C:\\Users\\yid\\TemporaryResearchDataStorage\\doe24_tspan_200\\StabilityAnalysisResults\\StabilityRegionsNutrient.csv",DataFrame)
df_microbial = CSV.read("C:\\Users\\yid\\TemporaryResearchDataStorage\\doe24_tspan_200\\StabilityAnalysisResults\\StabilityRegionsMicrobial.csv",DataFrame)
count_row = nrow(df_nutrient);
# count_feature = ncol(df_all) - 2
df1 = Matrix(df_microbial);
df2 = Matrix(df_nutrient);
A1 = df1[:,1:2];
Microbial_Region = df1[:,3];
A2 = df2[:,1:2];
Nutrient_Region = df2[:,3];
# println(Region)
count_total = 0;
count_0_region = 0;
count_1_region = 0;
count_2_region = 0;
New_data_0 = zeros(1,2)
New_data_1 = zeros(1,2)
New_data_2 = zeros(1,2)

for i in eachindex(Microbial_Region)
    count_total += 1
    # Go through all test data
    if Microbial_Region[i] == 0
        # println("Find parallel")
        count_0_region += 1
        # println(A[i,:])
        if count_0_region == 1
            New_data_0[:] = transpose(A1[i,:])
        else 
            New_data_0 = vcat(New_data_0, transpose(A1[i,:]))
        end

    elseif Microbial_Region[i] == 1
        count_1_region += 1
        if count_1_region == 1
            New_data_1[:] = transpose(A1[i,:])
        else 
            New_data_1 = vcat(New_data_1, transpose(A1[i,:]))
        end

    else Microbial_Region[i] == 2
        count_2_region += 1
        if count_2_region == 1
            New_data_2[:] = transpose(A1[i,:])
        else 
            New_data_2 = vcat(New_data_2, transpose(A1[i,:]))
        end
    end
end

scatter(New_data_0[:,1], New_data_0[:,2], xlabel = "C(Sucrose)", ylabel = "C(Ammonia)", label = "Exponential growth", markersize = 4, size = [1000,1000], markercolor = "Deep Sky Blue", markerstrokewidth = 0, legendfontsize=18, xtickfont = 18, ytickfont = 18, xguidefont = 18, yguidefont = 18)
scatter!(New_data_1[:,1], New_data_1[:,2], label = "With A Bump", markersize = 4,  markercolor = "Light Pink", markerstrokewidth = 0)
scatter!(New_data_2[:,1], New_data_2[:,2], label = "Without Bumps", markersize = 4, markercolor = "Dark Magenta", markerstrokewidth = 0)

Plots.savefig("C:\\Users\\yid\\TemporaryResearchDataStorage\\doe24_tspan_200\\StabilityAnalysisResults\\StabilityRegionsMicrobial.pdf")

count_total = 0;
count_0_region = 0;
count_1_region = 0;
count_2_region = 0;
New_data_0 = zeros(1,2)
New_data_1 = zeros(1,2)
New_data_2 = zeros(1,2)
for i in eachindex(Nutrient_Region)
    count_total += 1
    # Go through all test data
    if Nutrient_Region[i] == 0
        # println("Find parallel")
        count_0_region += 1
        # println(A[i,:])
        if count_0_region == 1
            New_data_0[:] = transpose(A2[i,:])
        else 
            New_data_0 = vcat(New_data_0, transpose(A2[i,:]))
        end
    
    # Monod L-P model
    elseif Nutrient_Region[i] == 2
        # println("Find parallel")
        count_2_region += 1
        # println(A[i,:])
        if count_2_region == 1
            New_data_2[:] = transpose(A2[i,:])
        else 
            New_data_2 = vcat(New_data_2, transpose(A2[i,:]))
        end

    else 
        count_1_region += 1
        if count_1_region == 1
            New_data_1[:] = transpose(A2[i,:])
        else 
            New_data_1 = vcat(New_data_1, transpose(A2[i,:]))
        end
    end
end
println(count_1_region)
scatter(New_data_0[:,1], New_data_0[:,2], xlabel = "C(Sucrose)", ylabel = "C(Ammonia)", label = "C Stable State", markersize = 4, size = [1000,1000], markercolor = "Deep Sky Blue", markerstrokewidth = 0, legendfontsize=18, xtickfont = 18, ytickfont = 18, xguidefont = 18, yguidefont = 18, legend=:topright)
scatter!(New_data_1[:,1], New_data_1[:,2], label = "C Not Touch Stable State", markersize = 4,  markercolor = "Light Pink", markerstrokewidth = 0)
scatter!(New_data_2[:,1], New_data_2[:,2], label = "N decreasing", markersize = 2, markercolor = "Dark Magenta", markerstrokewidth = 0)

Plots.savefig("C:\\Users\\yid\\TemporaryResearchDataStorage\\doe24_tspan_200\\StabilityAnalysisResults\\StabilityRegionsNutrient.pdf")

# using PlotlyJS
# plot([
#     scatter(x = New_data_0[:,1], y = New_data_0[:,2], fill = "tozeroy", mode = "none"),
#     scatter(x = New_data_1[:,1], y = New_data_1[:,2], fill = "tonexty", mode = "none"),
#     scatter(x = New_data_2[:,1], y = New_data_2[:,2], fill = "tonexty", mode = "none")
# ])
# plot([
#     scatter(x = New_data_0[:,1], y = New_data_0[:,2]),
#     scatter(x = New_data_1[:,1], y = New_data_1[:,2]),
#     scatter(x = New_data_2[:,1], y = New_data_2[:,2])
# ])

# plot([
#     scatter(x=1:4, y=[0, 2, 3, 5], fill="tozeroy"), # fill down to xaxis
#     scatter(x=1:4, y=[3, 5, 1, 7], fill="tonexty"), # fill to trace0 y
# ])
using VegaLite
df_all |>
stack |>

using CategoricalArrays
cv = categorical(df_all[!,"Stability"])
df_all[!,"RegionCategories"] = cv
@vlplot(:area, x=:Sucrose, y={:Ammonia, stack=:zero}, color = "RegionCategories:n")
df_all[!,:Stability] = parse.(Symbol, df_all[!, :Stability])
