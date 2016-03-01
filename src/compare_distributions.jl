immutable compare_distributions
    data::AbstractArray
    log_likehoods_ratio::AbstractArray
    sig_level::Float64
    xmin::Number
    V_test_stat::Float64
    V_p_val::Float64
    V_preff_distr::Int64
    C_b::Int64
    C_p_val::Float64
    C_preff_distr::Int64
end

function show(io::IO, x::compare_distributions)
    println(io, "Data:  $(x.data)")
    println(io, "Log likehood ratio: $(x.log_likehoods_ratio)")
    println(io, "Significance level: $(x.sig_level)")
    println(io, "Xmin: $(x.xmin)")
    println(io, "Vuong - test statictic: $(x.V_test_stat)")
    println(io, "Vuong - p-value: $(x.V_p_val)")
    println(io, "Vuong - preffered distribution: $(x.V_preff_distr)")
    println(io, "Clarke test - number of possitive values in log likehood ratio: $(x.C_b)")
    println(io, "Clarke test - p-value: $(x.C_p_val)")
    println(io, "Clarke test - preffered distribution: $(x.C_preff_distr)")
end



function compare_distributions(d1::con_powerlaw, d2::DataType, data::AbstractArray; sig_level = 0.05)
    if !(d2 <: ContinuousUnivariateDistribution)
        println("Both distributions should be continuos.")
        return Union{}
    end
    xmin = d1.θ
    data = sort(data)
    data = data[findfirst(x -> x >= xmin, data): end]
    d2 = fit(d2, data)
    
    _compare_distributions(d1,d2,data,xmin,sig_level)

end

function compare_distributions(d1::dis_powerlaw, d2::DataType, data::AbstractArray; sig_level = 0.05)
    if !(d2 <: DiscreteUnivariateDistribution)
        println("Both distributions should be discrete.")
        return Union{}
    end
    xmin = d1.θ
    data = sort(data)
    data = data[findfirst(x -> x >= xmin, data): end]
    d2 = fit(d2, data)
    
    _compare_distributions(d1,d2,data,xmin,sig_level)
end

function compare_distributions(d1::DataType, d2::DataType, data::AbstractArray, xmin::Number = 0; sig_level = 0.05)
    if !(((d1 <: ContinuousUnivariateDistribution) && (d2 <: ContinuousUnivariateDistribution)) || ((d1 <: DiscreteUnivariateDistribution) && (d2 <: DiscreteUnivariateDistribution)))
        println("Both distributions should be of same type (continuous or discrete).")
        return Union{}
    end

    if (xmin == 0)
        if (d1 <: DiscreteUnivariateDistribution) xmin = 1 end
        xmin = minimum(data)
    else
        data = sort(data)
        data = data[findfirst(x -> x >= xmin, data): end]
    end
    d1 = fit(d1, data)
    d2 = fit(d2, data)

    _compare_distributions(d1,d2,data,xmin,sig_level)
end

function compare_distributions(d1::DiscreteUnivariateDistribution, d2::DiscreteUnivariateDistribution, data::AbstractArray, xmin::Number = 1; sig_level = 0.05)
   _compare_distributions(d1,d2,data,xmin,sig_level) 
end

function compare_distributions(d1::ContinuousUnivariateDistribution, d2::ContinuousUnivariateDistribution, data::AbstractArray, xmin::Number = 0; sig_level = 0.05)
   _compare_distributions(d1,d2,data,xmin,sig_level) 
end

function _compare_distributions(d1::UnivariateDistribution, d2::UnivariateDistribution, data::AbstractArray, xmin::Number, sig_level)
    data = sort(data)
    xmin == 0 ? xmin = minimum(data) : data = data[findfirst(x -> x >= xmin, data): end]
    # Vuong's test
    log_likehoods_ratio = logpdf(d1,data) - logpdf(d2,data)
    n = length(log_likehoods_ratio)
    m = mean(log_likehoods_ratio)
    standard_deviation = std(log_likehoods_ratio)
    test_stat = sqrt(n) * m / standard_deviation
    v_p_val = cdf(Normal(), test_stat)
    v_preff_distr = 0
    if (test_stat > cdf(Normal(), 1 - sig_level/2))
        v_preff_distr = 1
    elseif (test_stat < -cdf(Normal(), 1 - sig_level/2))
        v_preff_distr = 2
    end
    
    # Clarke's test
    # b = number of possitive values in log_likehoods_ratio
    b = sum(x->x > 0 ? 1 : 0, log_likehoods_ratio)
    preff_distr = 0
    if (b >= n/2)
        pval = 2 * (1 - cdf(Binomial(n,0.5),b-1))
        if (pval <= sig_level) preff_distr = 1 end
    else
        pval = 2 * (cdf(Binomial(n,0.5),b))
        if (pval <= sig_level) preff_distr = 2 end
    end

    return compare_distributions(data,log_likehoods_ratio,sig_level,xmin,test_stat,v_p_val,v_preff_distr,b,pval,preff_distr)
end