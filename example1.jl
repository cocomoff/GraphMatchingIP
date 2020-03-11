# -*- coding: utf-8 -*-

using Cbc
using JuMP


function disp_mat(W)
    display(W); println()
end


# Synthetic graph weight
const L = 5 # bipartite 1 (左側のノード数)
const R = 5 # bipartite 2 (右側のノード数)
const W = round.(rand(Float32, L, R), digits=3) # ランダムな重み (wrd)


# スキップ処理 (適当に辺を消す．必要に応じて)
const is_skip = false
const θ = 0.2
if is_skip
    W[W .<= θ] .= 0
end
# disp_mat(W)

# parameters and clusters
εp, εl = 0.3, 0.5
K = 2 # G is partitioned into K sub-graphs


# Prepare link list
Nodes = collect(1:L+R)
Links = [(l, r + L) for l in 1:L, r in 1:R if W[l, r] > 0]
numV = L + R # number of nodes
numE = length(Links)
# @show Nodes
# @show Links


# solve IP
m = Model(with_optimizer(Cbc.Optimizer, LogLevel=0))

# variables (5f)
@variable(m, f[n=1:numV, k=1:K; n in Nodes], Bin) # node n is assigned to sub-grpah k
@variable(m, h[l=1:numE, k=1:K], Bin) # edge l is assigned to sub-graph k

# 目的関数 (5a)
obj = sum(W[Links[e][1], Links[e][2] - L] * sum(h[e, k] for k in 1:K) for e in 1:numE)
@objective(m, Max, obj)

# 制約 (5b) どこかに割り当てられる
for n in 1:numV
    @constraint(m, sum(f[n, k] for k in 1:K) == 1)
end

# 制約 (5c) 割当の巻き込み性質 (推移性)
for k in 1:K
    for e in 1:numE
        (p, pp) = Links[e]
        @constraint(m, f[p, k] + f[pp, k] >= 2 * h[e, k])
    end
end

# 制約 (5d) 割当のε性質その1
for k in 1:K
    @constraint(m, sum(h[e, k] for e in 1:numE) ≤ (1 + εl) * numE / K)
end

# 制約 (5e) 割当のε性質その2
for k in 1:K
    @constraint(m, sum(f[n, k] for n in 1:numV) ≤ (1 + εp) * numV / K)
end

# println(m)

# # 解く
optimize!(m)
obj_val = objective_value(m)
println(obj_val)

"""
@variable(m, f[n=1:numV, k=1:K; n in Nodes], Bin) # node n is assigned to sub-grpah k
@variable(m, h[l=1:numE, k=1:K], Bin) # edge l is assigned to sub-graph k
"""

# get assignment
for k in 1:K
    groupV = Set{Int}()
    groupE = Set{Tuple{Int, Int}}()
    for n in 1:numV
        (value(f[n, k]) > 0) && push!(groupV, n)
    end
    for l in 1:numE
        (value(h[l, k]) > 0) && push!(groupE, Links[l])
    end
    println("$k $groupV")
    for (l, r) in groupE
        println("  ($l, $r)")
    end
end