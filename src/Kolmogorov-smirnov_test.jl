function Kolmogorov_smirnov_test(dat::AbstractArray,d::ContinuousUnivariateDistribution,xmin::Number,xmax::Int64 = round(Int,1e5))
  data = sort(dat)
  n = float(length(data))
  max_indx = findlast(x -> x <= xmax, data)
  min_indx = findfirst(x -> x >= xmin, data)
  data = data[min_indx:max_indx]
  
  act_cdf = collect(0:length(data)-1) / n
  thr_cdf = cdf(d,float(data))
  D = maximum(abs(act_cdf-thr_cdf))
  return D
end

function Kolmogorov_smirnov_test(dat::AbstractArray,d::DiscreteUnivariateDistribution,xmin::Int64,xmax::Int64 = round(Int,1e5))
  alpha, xmin = params(d)
  data = round(Int,sort(dat))
  n = float(length(data))
  max_indx = findlast(x -> x <= xmax, data)
  min_indx = findfirst(x -> x >= xmin, data)
  data = data[min_indx:max_indx]
  
  thr_cdf = cdf(d,collect(xmin:data[end])+1)
  occurence = counts(data)[Int(xmin) - (data[1] - 1):end]
  act_cdf = occurence/sum(occurence)
  act_cdf = cumsum(act_cdf)
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
