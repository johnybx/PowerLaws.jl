module PowerLaws
  using StatsBase
  using Distributions
  using Optim
  using Compat
  import Distributions:rand,shape,pdf,ccdf,scale,params,cdf,cquantile,quantile,fit_mle,mean,median,var,skewness,mode,kurtosis,logccdf,logpdf,logcdf,entropy,@distr_support
  import Base:minimum,maximum,show

  export con_powerlaw,dis_powerlaw,Kolmogorov_smirnov_test,estimate_xmin,bootstrap,bootstrap_p,compare_distributions

  include("helper_macro.jl")
  include("powerlaw_discrete.jl")
  include("powerlaw_continuous.jl")
  include("estimate_xmin.jl")
  include("Kolmogorov-smirnov_test.jl")
  include("bootstrap.jl")
  include("compare_distributions.jl")

end



using Distributions
import PowerLaws:bootstrap,bootstrap_p,con_powerlaw,dis_powerlaw

#this parallel versions of bootstrap,bootstrap_p cannot be defined in module powerlaw at this time because
#of this issue #https://github.com/JuliaLang/julia/issues/13649

function bootstrap(data::AbstractArray,d::UnivariateDistribution, processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  α,θ = params(d)

  n = length(data)

  seed == 0 ? srand() : srand(seed)
  prcs = addprocs(processes)
  rmref = RemoteRef(1)
  put!(rmref,(d,n,xmins,xmax,data))
  helper_arr = Array(Int8,no_of_sims)
  fill!(helper_arr,0)
  @everywhere using PowerLaws
  @everywhere using Distributions
  helper_dict = Dict{Int64,RemoteRef}()
  for prc in workers() # copy data to all processes
    helper_dict[prc] = @spawnat prc begin
      global d,n,xmins,xmax,data
      d,n,xmins,xmax,data = fetch(rmref)
    end
  end
  for prc in workers() # wait till all data is coppied
    fetch(helper_dict[prc])
  end
  helper_dict = 0
  @everywhere function helper_func(args)
    sim_data = sample(data, n, replace=true)
    estimate_xmin(sim_data,typeof(d),xmins = xmins,xmax = xmax)
  end
  statistic = pmap(helper_func,helper_arr)
  rmprocs(prcs)
  take!(rmref)
  return statistic
end

function bootstrap(data::AbstractArray,distribution::Type{con_powerlaw},processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  d,ks = estimate_xmin(data,distribution,xmins = xmins,xmax =xmax)
  bootstrap(data,d,processes,no_of_sims = no_of_sims,xmins = xmins,xmax = xmax,seed =seed)
end

function bootstrap(data::AbstractArray,distribution::Type{dis_powerlaw},processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  d,ks = estimate_xmin(data,distribution,xmins = xmins,xmax =xmax)
  bootstrap(data,d,processes,no_of_sims = no_of_sims,xmins = xmins,xmax = xmax,seed =seed)
end

function bootstrap_p(data::AbstractArray,d::UnivariateDistribution, processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  sort_data = sort(data)
  α,θ = params(d)
  n = length(sort_data)
  tail_indx = findfirst(sort_data,θ)
  tail_p = length(sort_data[tail_indx:end])/n
  KS_stat = Kolmogorov_smirnov_test(sort_data[tail_indx:end],d)

  seed == 0 ? srand() : srand(seed)
  prcs = addprocs(processes)
  rmref = RemoteRef(1)
  put!(rmref,(d,n,tail_p,tail_indx,xmins,xmax,sort_data))
  helper_arr = Array(Int8,no_of_sims)
  fill!(helper_arr,0)
  @everywhere using PowerLaws
  @everywhere using Distributions

  helper_dict = Dict{Int64,RemoteRef}()
  for prc in workers() # copy data to all processes
    helper_dict[prc] = @spawnat prc begin
      global d,n,tail_p,tail_indx,xmins,xmax,sort_data
      d,n,tail_p,tail_indx,xmins,xmax,sort_data = fetch(rmref)
    end
  end
  for prc in workers() # wait till all data is coppied
    fetch(helper_dict[prc])
  end
  helper_dict = 0

  @everywhere helper_func(args) = begin
    n1 = sum(map(x-> x>tail_p,rand(n)))
    n2 = n - n1
    sim_data = Array(Float64,0)
    append!(sim_data,sample(sort_data[1:tail_indx-1],n1,replace = true))
    append!(sim_data,rand(d,n2))
    estimate_xmin(sim_data,typeof(d),xmins = xmins,xmax = xmax)
  end

  statistic= pmap(helper_func,helper_arr)

  rmprocs(prcs)
  take!(rmref)
  P = 0
  for i=1:no_of_sims
    if (KS_stat <= statistic[i][2])
      P +=1
    end
  end

  return statistic,(P/no_of_sims)
end

function bootstrap_p(data::AbstractArray,distribution::Type{con_powerlaw},processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  d,ks = estimate_xmin(data,distribution,xmins = xmins,xmax =xmax)
  bootstrap_p(data,d,processes,no_of_sims = no_of_sims,xmins = xmins,xmax = xmax,seed =seed)
end

function bootstrap_p(data::AbstractArray,distribution::Type{dis_powerlaw},processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
  d,ks = estimate_xmin(data,distribution,xmins = xmins,xmax =xmax)
  bootstrap_p(data,d,processes,no_of_sims = no_of_sims,xmins = xmins,xmax = xmax,seed =seed)
end
