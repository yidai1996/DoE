# Kinetic modeling and stability analysis of mutual cross feeding
This project develops a modeling framework for the case of a two-strain mutual cross-feeding system where the consumption of a substrate by one of the strains produces metabolites that can act as resources for the other. 

---

## Quick start

```bash
git clone https://github.com/yidai1996/DoE.git
cd DoE
# Activate the Julia environment and install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run the main script
julia --project=. Models/Monoculture_producer_TriCulture_Chemostat_Monod.jl
