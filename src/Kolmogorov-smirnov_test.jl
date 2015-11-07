function Kolmogorov_smirnov_test(dat::AbstractArray,d::con_powerlaw,xmax::Int64 = round(Int,1e5))
  alpha, xmin = params(d)
  data = sort(dat)
  n = float(length(data))
  if (data[end] > xmax)
    for i=1:length(data)
      if (data[end - i]<=xmax)
        data = data[1:end - i]
        break
      end
    end
  end
  if (length(data) < 2)
    println("Data must be at least of lenght 2!")
    return Inf
  end
  act_cdf = [0:length(data)-1] / n
  thr_cdf = cdf(d,float(data))
  D = maximum(abs(act_cdf-thr_cdf))
  return D
end

function Kolmogorov_smirnov_test(dat::AbstractArray,d::dis_powerlaw,xmax::Int64 = round(Int,1e5))
  alpha, xmin = params(d)
  data = sort(dat)
  n = float(length(data))
  if (data[end] > xmax)
    for i=1:length(data)
      if (data[end - i]<=xmax)
        data = data[1:end - i]
        break
      end
    end
  end
  if (length(data) < 2)
    println("Data must be at least of lenght 2!")
    return Inf
  end

  unique_data = sort(unique(data))
  element = data[1]
  indx = 0
  arr = Array(Float64,0)
  for i=1:length(data)
    if element == data[i]
      push!(arr,indx)
    else
      indx=i-1
      element = data[i]
      push!(arr,indx)
    end
  end
  arr = sort(unique(arr))
  act_cdf = arr/n
  if (1-zeta(alpha,xmin) == 1)
    thr_cdf = ones(Float64,length(unique_data))
  else
    thr_cdf = cdf(d,unique_data)
  end

  D = maximum(abs(act_cdf-thr_cdf))
  return D
end

#helper function
function create_histogram(x::AbstractArray)
  h = zeros(typeof(x[1]),0)
  max = 0
  d = 0
  for i = 1:length(x)
    if max < x[i]
      for j=1:(x[i]-max)
        push!(h,0)
      end
      h[x[i]] =1
      max = x[i]
    else
      h[x[i]] +=1
    end
  end
  return h
end
