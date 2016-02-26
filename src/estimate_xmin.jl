function init_xmins(data::AbstractArray,xmins::AbstractArray,xmax::Int64)
  sorted_data = sort(data)
  bins_data = array_bins(sorted_data)
  if (xmins == [])
    xmins = sort(unique(sorted_data))
    if (xmins[end] > xmax)
      max_indx = findlast(x -> x <= xmax, xmins)
      xmins = xmins[1:max_indx]
    end
    xmins = xmins[1:end]
    if xmins[1] <= 0
      println("removing elements smaller than 0")
      xmins = xmins[findfirst(x -> x > 0,xmins), end]
    end
  else
    real_xmins = Array(Bool,length(xmins))
    for i=1:length(xmins)
      if (haskey(bins_data,xmins[i]))
        real_xmins[i] = true
      else
        real_xmins[i] = false
      end
    end
    xmins = xmins[real_xmins]
  end
  return sorted_data,bins_data,xmins
end

function _estimate_xmin(sorted_data::AbstractArray,bins_data::Dict{Float64,Int64},distribution::DataType,xmins::AbstractArray,xmax::Int64)
  if (length(xmins) == 0)
    println("No xmins")
    return Union{}
  end

  min_dist = Inf
  best_fit = Union{}

  for xmin in xmins
    fit_data = sorted_data[bins_data[xmin]:end]
    f = fit(distribution,fit_data)
    
    negloglike(alpha) = begin
      d = distribution(alpha[1],f.θ)
      r = -sum(logpdf(d,fit_data))
      if(Inf == r || -Inf == r)
        r = 1e12
      end
      return r
    end
    try
      opt_alfa = fminbox(DifferentiableFunction(negloglike),[f.α],[1.0],[Inf])
      f = distribution(opt_alfa.minimum[1],f.θ)
    catch
      #if fminbox throws error it means that function cannot be optimized
    end
    if (f.α == Inf)
      continue
    end
    d = Kolmogorov_smirnov_test(fit_data,f,xmin,xmax)

    if ((min_dist > d))
      best_fit = f
      min_dist = d
    end
  end
  return best_fit,min_dist
end

function estimate_xmin(data::AbstractArray,distribution::Type{con_powerlaw};xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5))

  sorted_data,bins_data,xmins = init_xmins(data,xmins,xmax)

  _estimate_xmin(sorted_data,bins_data,distribution,xmins,xmax)
end

function estimate_xmin(data::AbstractArray,distribution::Type{dis_powerlaw};xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5))
  if (data != floor(data))
    println("Data should be discreate. Use round or floor function.")
    return Union{}
  end

  sorted_data,bins_data,xmins = init_xmins(data,xmins,xmax)
  if xmins[1] < 1
    xmins = xmins[findfirst(x -> x >= 1), end]
  end

  _estimate_xmin(sorted_data,bins_data,distribution,xmins,xmax)
end

function estimate_xmin(data::AbstractArray,distribution::DataType;xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5))
  min_dist = Inf
  best_fit = Union{}
  if ((distribution <: DiscreteUnivariateDistribution)&&(data != floor(data)))
    println("Data should be discreate. Use round or floor function.")
    return Union{}
  end

  sorted_data,bins_data,xmins = init_xmins(data,xmins,xmax)
  if ((distribution <: DiscreteUnivariateDistribution)&&(xmins[1] < 1))
    xmins = xmins[findfirst(x -> x >= 1), end]
  end

  if (length(xmins) == 0)
    println("No xmins")
    return Union{}
  end

  for xmin in xmins
    fit_data = sorted_data[bins_data[xmin]:end]
    f = fit(distribution,fit_data)
    d = Kolmogorov_smirnov_test(fit_data,f,xmin,xmax)

    if ((min_dist > d))
      best_fit = f
      min_dist = d
    end
  end
  return best_fit,min_dist
end

#helper function  return dict{Int,Int} where key is number and value is where in array it starts
#input array must be sorted
function array_bins(arr::AbstractArray)
  num = arr[1]
  bins = Dict{Float64,Int64}()
  bins[num] = 1
  for i=1: length(arr)
    if num != arr[i]
      num = arr[i]
      bins[num] = i
    end
  end
  return bins
end
