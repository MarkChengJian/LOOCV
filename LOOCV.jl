using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Distributions")
Pkg.add("Plots")
Pkg.add("LinearAlgebra")
Pkg.add("Statistics")
Pkg.add("Random")
Pkg.add("DelimitedFiles")

using Distributions
using LinearAlgebra
using Random
using CSV
using DataFrames


# Efficient way for LOOCV
function loiEBVEff2(i,Vi,C,y)
    r     = 1/Vi[i,i]
    q     = -Vi[i,:] *r
    yi    = copy(y)
    yi[i] = 0.0
    Ci    = C[i,:]
    Ci[i] = 0.0
    (Ci'*Vi - (Ci'q) * q'/r)*yi
end



function LOOCV(genotypes, phenotypes, t, Ve, Vg, sum2pq, output)
     res = phenotypes[:,t]
     sel = .!(ismissing.(res))
     ZFull = Matrix{Float64}(I,nrow(phenotypes),nrow(phenotypes))
     Z = ZFull[sel,:]
    
    varAlpha = Vg/sum2pq
    G = Z*M2*M2'*Z'*varAlpha
    V = G + I*Ve
    
    Vi = inv(V)
    C = G
    n = size(V,1)
    phenotypes=dropmissing(phenotypes, t)
    y  = phenotypes[!,t]
    
    EBV2=[loiEBVEff2(i,Vi,C,y) for i=1:n]
    print(cor(EBV2, y))
    
    ID= phenotypes[:,"germplasm_id"]
    YR= phenotypes[:,"PS3_YEAR"]
    result=DataFrame([ID YR EBV2 y],:auto)
    CSV.write(output, result)
    
end


