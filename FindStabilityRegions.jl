using XLSX, DataFrames, CSV

# Distinguish stability by checking the last two points in all growth curves
out_dir = "C:\\Users\\yid\\TemporaryResearchDataStorage\\doe3_tspan_200\\excel\\"
index_files = readdir(out_dir);
# build a dataframe for C/N and stability (growth curve regions)
df = DataFrame([Int64[], Int64[], Int64[]], ["C(Sucrose)","Ammonia","Stability"])
for i in eachindex(index_files)
    # search all excel files
    # find the last two points of biomass concentration
    ex = XLSX.readxlsx(out_dir * index_files[i])
    sh = ex["Sheet1"]
    num_rows, num_cols = size(sh[:])
    Av1 = sh["B" * string(num_rows-1)]
    Av2 = sh["B" * string(num_rows)]
    Se1 = sh["C" * string(num_rows-1)]
    Se2 = sh["B" * string(num_rows)]
    println("Av=",[Av1 Av2])
    println("Se=",[Se1 Se2])
    # Check the stability
    if num_rows<=3
        # the stable region in linear stability analysis
        push!(df,[sh["D2"], sh["E2"], 0])
    elseif (Av1-Av2 <= 10^(-5)) & (Se2-Se1 <= 10^(-5))
        # going to exponentially growth
        push!(df,[sh["D2"], sh["E2"], 1])
    else push!(df,[sh["D2"], sh["E2"], 2])
    end
    println(i)
end
println(df)
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\doe3_tspan_200\\StabilityAnalysisResults\\StabilityRegions.csv",df)