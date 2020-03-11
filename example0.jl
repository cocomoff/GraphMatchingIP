# -*- coding: utf-8 -*-

using Cbc
using JuMP


function disp_mat(W)
    display(W); println()
end


# Synthetic graph weight
const L = 10 # bipartite 1 (左側のノード数)
const R = 10 # bipartite 2 (右側のノード数)
const W = rand(Float32, L, R) # ランダムな重み (wrd)

# スキップ処理 (適当に辺を消す．必要に応じて)
const is_skip = false
const θ = 0.2
if is_skip
    W[W .<= θ] .= 0
end
# disp_mat(W)

# solve IP
m = Model(with_optimizer(Cbc.Optimizer, LogLevel=0))

@variable(m, u[l=1:L, r=1:R], Bin)  # 変数 (4d)
for l in 1:L
    @constraint(m, sum(u[l, r] for r in 1:R) ≤ 1)  # 制約 (4c)
end
for r in 1:R
    @constraint(m, sum(u[l, r] for l in 1:L) ≤ 1)  # 制約 (4b)
end

@objective(m, Max, sum(W[l, r] * u[l, r] for l in 1:L, r in 1:R if W[l, r] > 0)) # 目的関数 (4a)

# println(m)

# 解く
optimize!(m)
obj_val = objective_value(m)
println(obj_val)

for l in 1:L, r in 1:R
    if value(u[l, r]) > 0
        println("selected (l,r)=($l,$r)")
    end
end