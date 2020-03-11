# -*- coding: utf-8 -*-

using Cbc
using JuMP


# Detour Weights (TODO: need to compute here properly)
const R = 5
const L = 5
const N = R + L # number of trips
const C = round.(rand(Float32, N, N), digits=3)

# parameters and clusters
εp, εl = 0.3, 0.5
ε = 0.3
K = 2 # G is partitioned into K sub-graphs


# Model

m = Model(with_optimizer(Cbc.Optimizer, LogLevel=0))

# decision variables (8f, 8g)
@variable(m, q[n=1:N, k=1:N], Bin)
@variable(m, ζ[k=1:N], Bin)

# objective (8a)
@objective(m, Min, sum(C[n, k] * q[n, k] for n in 1:N, k in 1:N))

# 制約 (8b)
@constraint(m, sum(ζ[k] for k in 1:N) == K)

# 制約 (8c)
for n in 1:N
    @constraint(m, sum(q[n, k] for k in 1:N) == 1)
end

# 制約 (8d)
for k in 1:N
    @constraint(m, sum(q[n, k] for n in 1:R) ≤ (1 + ε) * R / K * ζ[k])
end

# 制約 (8e)
for k in 1:N
    @constraint(m, sum(q[n, k] for n in R+1:N) ≤ (1 + ε) * L/  K * ζ[k])
end


# # 解く
optimize!(m)
obj_val = objective_value(m)
println(obj_val)

"""
@variable(m, q[n=1:N, k=1:N], Bin)
@variable(m, ζ[k=1:N], Bin)
"""

# get representative
for k in 1:N
    (value(ζ[k]) > 0) && println("$k $(value(ζ[k]))")
end

# get values in q
for n in 1:N, k in 1:N
    (value(q[n, k]) > 0) && println("(n, k)=($n, $k)")
end