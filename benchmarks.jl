using DataFrames, StatFiles, Inequality, BenchmarkTools, StatsBase, Bootstrap, 
    Random, CSV, Tidier

"""
    preprocess_trial(t::BenchmarkTools.Trial, id::AbstractString)

Extract useful information from benchmark trials. 
"""
preprocess_trial(t::BenchmarkTools.Trial, id::AbstractString) =
    (id=id,
     median=median(t.times),
     allocations=t.allocs,
     memory_estimate=t.memory)


"""
    compute_groupby_estimate_with_boot(dataframe, groupvar, func, weights, m)

Compute an estimate with bootstrap and grouped by a variable.
"""
function compute_groupby_estimate_with_boot(dataframe, groupvar, func, weights, m)::DataFrame
    return combine(groupby(dataframe, Symbol(groupvar)),
    df -> combine(df, df_ -> begin
        return bootstrap(d -> func(d[:,:dhi], Weights(d[:,weights])), df_,  BasicSampling(m))
    end)[1,1].t0[1])
end 


"""
    recode_hhtype(hhtype, nhhmem17)

Recode the 'hhtype' variable into categories (1 to 5). 
"""
function recode_hhtype(hhtype, nhhmem17)::Int64
    if hhtype == 100
        return 1
    elseif hhtype == 210
        return 2
    elseif (hhtype == 220) && nhhmem17 > 1
        return 3
    elseif (hhtype == 230) && nhhmem17 > 1
        return 4
    elseif 230 < hhtype <= 930
        return 5
    elseif (hhtype == 220) && nhhmem17 == 0
        return 5
    elseif (hhtype == 230) && nhhmem17 == 0 
        return 5
    else 
        return 6
    end
end


# maximum number of seconds in each benchmark
max_seconds = 1200
bootstrap_n = 1000

Random.seed!(4)

df_raw = DataFrame(load("/home/data/it14ih.dta"))

# df_raw = df_raw[:,[:hid,:dhi,:hwgt]]

df = df_raw

df = df[:,[:hid,:dhi,:hwgt]]

df_0 = df

# Enlarge the DF (to 102,000 obs)
for i in range(1, 101)
    global df
    df = vcat(df, df_0)
end

dropmissing!(df, disallowmissing=true) # 100,062 remaining

# Change variable type so it can be stored
df.hid = convert(Vector{Int64}, df.hid)

df = subset(df, :dhi => ByRow(x -> x > 0 ) )


dhi = df.dhi
hwgt = df.hwgt

# Benchmarks



benchmarks = DataFrame()



# Simple estimates --------------------------------------------------------

println("Starting the computation of benchmarks for simple estimates")

# gini
benchmark_gini = @benchmark gini(dhi, hwgt)
push!(benchmarks, preprocess_trial(benchmark_gini, "gini"))


# atkinson with ϵ > 1
benchmark_atkinson_large_ϵ = @benchmark atkinson(dhi, hwgt, 1.2)
push!(benchmarks, preprocess_trial(benchmark_atkinson_large_ϵ, "benchmark_atkinson_large_ϵ"))


# atkinson with ϵ < 1
benchmark_atkinson_small_ϵ = @benchmark atkinson(dhi, hwgt, 0.8)
push!(benchmarks, preprocess_trial(benchmark_atkinson_small_ϵ, "benchmark_atkinson_small_ϵ"))


# Fgt
benchmark_fgt = @benchmark fgt(dhi, hwgt, 0.8, median(dhi, pweights(hwgt)) .* 0.6 )
push!(benchmarks, preprocess_trial(benchmark_fgt, "benchmark_fgt"))


# headcount
benchmark_headcount = @benchmark headcount(dhi, hwgt, median(dhi, pweights(hwgt)) .* 0.6)
push!(benchmarks, preprocess_trial(benchmark_headcount, "benchmark_headcount"))


# poverty gap
benchmark_poverty_gap = @benchmark poverty_gap(dhi, hwgt, median(dhi, pweights(hwgt)) .* 0.6)
push!(benchmarks, preprocess_trial(benchmark_poverty_gap, "benchmark_poverty_gap"))


# watts
benchmark_watts = @benchmark watts(dhi, hwgt, median(dhi, pweights(hwgt)) .* 0.6)
push!(benchmarks, preprocess_trial(benchmark_watts, "benchmark_watts"))


# theil
benchmark_theil = @benchmark theil(dhi, hwgt)
push!(benchmarks, preprocess_trial(benchmark_theil, "benchmark_theil"))


# # lorenz_curve
# benchmark_lorenz_curve = @benchmark lorenz_curve(dhi, hwgt)
# push!(benchmarks, preprocess_trial(benchmark_lorenz_curve, "benchmark_lorenz_curve"))


# mld
benchmark_mld = @benchmark mld(dhi, hwgt)
push!(benchmarks, preprocess_trial(benchmark_mld, "benchmark_mld"))


# Estimates with bootstrap ------------------------------------------------

println("Starting the computation of benchmarks for estimates with bootstrap")

# gini
benchmark_bootstrap_gini = @benchmark bootstrap(d -> gini(d[:,:dhi], d[:,:hwgt]), df, BasicSampling(bootstrap_n)) samples=100 seconds=max_seconds
push!(benchmarks, preprocess_trial(benchmark_bootstrap_gini, "benchmark_bootstrap_gini"))


# atkinson with ϵ > 1
benchmark_bootstrap__atkinson_large_ϵ = @benchmark bootstrap(d -> atkinson(d[:,:dhi], d[:,:hwgt], 1.2), df, BasicSampling(bootstrap_n)) samples=100 seconds=max_seconds
push!(benchmarks, preprocess_trial(benchmark_bootstrap__atkinson_large_ϵ, "benchmark_bootstrap__atkinson_large_ϵ"))


# atkinson with ϵ < 1
benchmark_bootstrap__atkinson_small_ϵ = @benchmark bootstrap(d -> atkinson(d[:,:dhi], d[:,:hwgt], 0.8), df, BasicSampling(bootstrap_n)) samples=100 seconds=max_seconds
push!(benchmarks, preprocess_trial(benchmark_bootstrap__atkinson_small_ϵ, "benchmark_bootstrap__atkinson_small_ϵ"))


# Fgt
benchmark_bootstrap_fgt = @benchmark bootstrap(d -> fgt(d[:,:dhi], d[:,:hwgt], 0.8, median(dhi, pweights(hwgt)) .* 0.6), df, BasicSampling(bootstrap_n)) samples=100 seconds=max_seconds
push!(benchmarks, preprocess_trial(benchmark_bootstrap_fgt, "benchmark_bootstrap_fgt"))


# headcount
benchmark_bootstrap_headcount = @benchmark bootstrap(d -> headcount(d[:,:dhi], d[:,:hwgt], median(d[:,:dhi], pweights(d[:,:hwgt])) .* 0.6), df, BasicSampling(bootstrap_n)) samples=100 seconds=max_seconds
push!(benchmarks, preprocess_trial(benchmark_bootstrap_headcount, "benchmark_bootstrap_headcount"))


# poverty gap
benchmark_bootstrap_gap = @benchmark bootstrap(d -> poverty_gap(d[:,:dhi], d[:,:hwgt], median(d[:,:dhi], pweights(d[:,:hwgt])) .* 0.6), df, BasicSampling(bootstrap_n)) samples=100 seconds=max_seconds
push!(benchmarks, preprocess_trial(benchmark_bootstrap_gap, "benchmark_bootstrap_gap"))


# watts
benchmark_bootstrap_watts = @benchmark bootstrap(d -> watts(d[:,:dhi], d[:,:hwgt], 1.2), df, BasicSampling(bootstrap_n)) samples=100 seconds=max_seconds
push!(benchmarks, preprocess_trial(benchmark_bootstrap_watts, "benchmark_bootstrap_watts"))


# theil
benchmark_bootstrap_theil = @benchmark bootstrap(d -> theil(d[:,:dhi], d[:,:hwgt]), df, BasicSampling(bootstrap_n)) samples=100 seconds=max_seconds
push!(benchmarks, preprocess_trial(benchmark_bootstrap_theil, "benchmark_bootstrap_theil"))


# mld
benchmark_bootstrap_mld = @benchmark bootstrap(d -> mld(d[:,:dhi], d[:,:hwgt]), df, BasicSampling(bootstrap_n)) samples=100 seconds=max_seconds
push!(benchmarks, preprocess_trial(benchmark_bootstrap_mld, "benchmark_bootstrap_mld"))


# By groups ---------------------------------------------------------------

# use hhtype to group by

df = df_raw[:,[:hid,:dhi,:hhtype,:nhhmem17,:hwgt]]

dropmissing!(df, disallowmissing=true) 


df_0 = df

# Enlarge the DF (to 102,000 obs)
for i in range(1, 101)
    global df
    df = vcat(df, df_0)
end

df = subset(df, :dhi => ByRow(x -> x > 0 ) )



# ageyoch_less18 = if_else(any(ageyoch < 18) | (unique(nhhmem17) > 0) , TRUE, FALSE, missing=FALSE)

# 100 = 1
# 210 = 2
# 220 and ageyoch_less18 = 3
# 230 and ageyoch_less18 = 4
# between(220, 930) = 5
# hhtype 220 and not ageyoch_less18 = 5
# hhtype 230 and not ageyoch_less18 = 5

transform!(df, [:hhtype,:nhhmem17] => ByRow((A, B) -> recode_hhtype(A, B)) => :htype)

# @chain df_raw begin
#     @mutate(htype = recode_hhtype(hhtype, nhhmem17))
# end

benchmark_bootstrap_grouped_mld = @benchmark compute_groupby_estimate_with_boot(df, "htype", mld, "hwgt", 1000) samples=100 seconds=60
push!(benchmarks, preprocess_trial(benchmark_bootstrap_grouped_mld, "benchmark_bootstrap_grouped_mld"))


# Output benchmarks -------------------------------------------------------

CSV.write("/home/outputs/benchmarks_julia.csv", benchmarks)