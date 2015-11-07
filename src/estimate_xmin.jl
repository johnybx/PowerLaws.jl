function estimate_xmin(data::AbstractArray,discrete::Bool = false;xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5))
  min_dist = Inf
  best_fit = None #Union{} for v0.4
  sorted_data = sort(data)
  bins_data = array_bins(sorted_data)
  if (xmins == [])
    xmins = sort(unique(sorted_data))
    if (xmins[end] > xmax)
      for i=1:length(xmins)
        if (xmins[end-i] <= xmax)
          xmins = xmins[1:end-i]
          break
        end
      end
    end
    xmins = xmins[1:end-2]
    if xmins[1] < 1
      println("removing elements smaller than 1")
      while xmins[1] < 1
        shift!(xmins)
      end
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
  if (length(xmins) == 0)
    println("No xmins")
    return None #Union{}
  end
  for xmin in xmins
    fit_data = sorted_data[bins_data[xmin]:end]
    if discrete
      f = fit(dis_powerlaw,fit_data)

      negloglike(alpha) = begin
        d = dis_powerlaw(alpha[1],f.θ)
        r = -sum(logpdf(d,fit_data))
        if(Inf == r || -Inf == r)
          r = 1e12
        end
        return r
      end
      try
        opt_alfa = fminbox(DifferentiableFunction(negloglike),[f.α],[1.0],[Inf])
        f = dis_powerlaw(opt_alfa.minimum[1],f.θ)
      catch
        #if fminbox throws error it means that function cannot be optimized
      end

    else
      f = fit(con_powerlaw,fit_data)
    end
    d = Kolmogorov_smirnov_test(fit_data,f,xmax)

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
