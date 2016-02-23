function bootstrap(data::AbstractArray,d::UnivariateDistribution;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  n = length(data)

  seed == 0 ? srand() : srand(seed)
  statistic = Array(Tuple{typeof(d),Float64},no_of_sims)
  for i=1:no_of_sims
    sim_data = sample(data, n, replace=true)
    statistic[i] =estimate_xmin(sim_data,typeof(d),xmins = xmins,xmax = xmax)
  end
  return statistic
end

function bootstrap(data::AbstractArray,distribution::Type{con_powerlaw};no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  d,ks = estimate_xmin(data,distribution,xmins = xmins,xmax =xmax)
  bootstrap(data,d,no_of_sims = no_of_sims,xmins = xmins,xmax = xmax,seed =seed)
end

function bootstrap(data::AbstractArray,distribution::Type{dis_powerlaw};no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  d,ks = estimate_xmin(data,distribution,xmins = xmins,xmax =xmax)
  bootstrap(data,d,no_of_sims = no_of_sims,xmins = xmins,xmax = xmax,seed =seed)
end

function bootstrap_p(data::AbstractArray,d::UnivariateDistribution;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  sort_data = sort(data)
  α,θ = params(d)
  n = length(sort_data)
  tail_indx = findfirst(sort_data,θ)
  tail_p = length(sort_data[tail_indx:end])/n
  KS_stat = Kolmogorov_smirnov_test(sort_data[tail_indx:end],d)

  P = 0
  statistic = Array(Tuple{typeof(d),Float64},no_of_sims)
  seed == 0 ? srand() : srand(seed)
  for i=1:no_of_sims
    n1 = sum(map(x-> x>tail_p,rand(n)))
    n2 = n - n1
    sim_data = Array(Float64,0)
    append!(sim_data,sample(sort_data[1:tail_indx-1],n1,replace = true))
    append!(sim_data,rand(d,n2))
    statistic[i] = estimate_xmin(sim_data,typeof(d),xmins = xmins,xmax = xmax)
    if (KS_stat <= statistic[i][2])
      P +=1
    end
  end
  return statistic,(P/no_of_sims)
end

function bootstrap_p(data::AbstractArray,distribution::Type{con_powerlaw};no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  d,ks = estimate_xmin(data,distribution,xmins = xmins,xmax =xmax)
  bootstrap_p(data,d,no_of_sims = no_of_sims,xmins = xmins,xmax = xmax,seed =seed)
end

function bootstrap_p(data::AbstractArray,distribution::Type{dis_powerlaw};no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  d,ks = estimate_xmin(data,distribution,xmins = xmins,xmax =xmax)
  bootstrap_p(data,d,no_of_sims = no_of_sims,xmins = xmins,xmax = xmax,seed =seed)
end
