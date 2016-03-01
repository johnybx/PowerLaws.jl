Documentation
========

### Type definition

This definitions are done according to [Distributions](http://distributionsjl.readthedocs.org/en/latest/) package template. Continuous powerlaw support all function like cdf,ccdf,pdf,median,var etc. Discrete powerlaw do not have implemented moments - median,variance etc.

#### con_powerlaw(alfa,xmin)
Continuios power law with probability density function ![equation](http://www.sciweavers.org/tex2img.php?eq=p%28x%29%20%3D%5Cbegin%7Bcases%7D%20%5Cfrac%7B%28%5Calpha%20%20-%201%29%7D%7B%5Ctheta%7D%20%2A%20%28%5Cfrac%7Bx%7D%7B%5Ctheta%7D%29%5E%7B-%5Calpha%7D%20%20%20%26%20x%20%20%5Cgeq%20%20%5Ctheta%20%5C%5C0%20%26%20x%20%20%3C%20%5Ctheta%20%5Cend%7Bcases%7D%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0%22%20align=%22center%22%20border=%220%22%20alt=%22p(x)%20=\begin{cases}%20\frac{(\alpha%20%20-%201)}{\theta}%20*%20(\frac{x}{\theta})^{-\alpha}%20%20%20&%20x%20%20\geq%20%20\theta%20\\0%20&%20x%20%20%3C%20\theta%20\end{cases}%20%22%20width=%22244%22%20height=%2250%22)

#### dis_powerlaw(alfa,xmin)
Discrete power law with probability mass function ![equation](http://www.sciweavers.org/tex2img.php?eq=p%28x%29%20%3D%5Cbegin%7Bcases%7D%20%5Cfrac%7Bx%5E%7B-%5Calpha%7D%7D%7Bzeta%28%20%5Calpha%20%2C%20%5Ctheta%20%29%7D%20%20%20%26%20x%20%20%5Cgeq%20%20%5Ctheta%20%5C%5C0%20%26%20x%20%20%3C%20%5Ctheta%20%5Cend%7Bcases%7D%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0%22%20align=%22center%22%20border=%220%22%20alt=%22p(x)%20=\begin{cases}%20\frac{x^{-\alpha}}{zeta(%20\alpha%20,%20\theta%20)}%20%20%20&%20x%20%20\geq%20%20\theta%20\\0%20&%20x%20%20%3C%20\theta%20\end{cases}%20%22%20width=%22201%22%20height=%2254%22)

#### immutable compare_distributions (data::AbstractArray, log_likehoods_ratio::AbstractArray, sig_level::Float64, xmin::Number, V_test_stat::Float64, V_p_val::Float64, V_preff_distr::Int64, C_b::Int64, C_p_val::Float64, C_preff_distr::Int64,)

Data type for distribution comparison.

**Parameters:**

* data: array which according which should be distributions compared
* log_likehoods_ratio: log likehood ratio of data
* sig_level: sigma level
* xmin: smallest element which was used for comparing distributions
* V_test_stat: Vuong test statistic
* V_p_val: p-value from Vuong test
* V_preff_distr: preffered distibution according to Vuong test
* C_b: Clarke test b value - sum of possitive values in log likehood ratio
* C_p_val: p-value from Clarke test
* C_preff_distr: preffered distibution according to Clarke test

### Functions

#### estimate_xmin(data::AbstractArray,distribution::Type{con_powerlaw},xmins::AbstractArray = [],xmax::Int64 = int(1e5))
#### estimate_xmin(data::AbstractArray,distribution::Type{dis_powerlaw},xmins::AbstractArray = [],xmax::Int64 = int(1e5))

Estimates best xmin and alfa for data set with respect to the [Kolmogorov smirnov test](https://www.encyclopediaofmath.org/index.php/Kolmogorov-Smirnov_test)

**Parameters:**

* data: array of data which should be fit to distribution
* distribution: this parameter sets distribution Type at this time only dis_powerlaw,con_powerlaw
* xmins: if not specified all unique values in data are taken as possible xmins, if specified then only values in array xmins are considered when finding best xmin
* xmax: maximum value considered in calculations,values above xmax are not considered in for example calculating Kolmogorov smirnov test.

return powerlaw distribution (continuos or discrete) with best xmin and alfa and Kolmogorov smirnov test which belongs to the distribution.

**Examples**

````jl
    estimate_xmin(collect(1:100),con_powerlaw)
    estimate_xmin(collect(1:100),dis_powerlaw)
    estimate_xmin(collect(1:100),con_powerlaw,xmins = [1,2,3,4,50,90])

````
    
#### Kolmogorov_smirnov_test(dat::AbstractArray,d::con_powerlaw,xmax::Int64 = int(1e5))
#### Kolmogorov_smirnov_test(dat::AbstractArray,d::dis_powerlaw,xmax::Int64 = int(1e5))

Calculate [Kolmogorov smirnov test](https://www.encyclopediaofmath.org/index.php/Kolmogorov-Smirnov_test) on given data and distribution.

**Parameters:**

* data: array of data which should be fit to distribution
* d: distribution, right now it support only con_powerlaw and dis_powerlaw
* xmax: maximum value considered in calculations,values above xmax are not considered.

return KS_statistic

**Examples**

````jl
    Kolmogorov_smirnov_test([1:100],con_powerlaw(2,2))
````

#### bootstrap(data::AbstractArray,d::UnivariateDistribution;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
To quantify the uncertainty in our estimate for xmin you can use bootstrap method. More information can be found in this document [Power-law distributions in empirical data](http://arxiv.org/pdf/0706.1062v2.pdf)

**Parameters:**

* data: array of data which should be fit to distribution
* d: distribution, right now it support only con_powerlaw and dis_powerlaw
* no_of_sims: number of simulations which should be performed
* xmins: if not specified all unique values in data are taken as possible xmins, if specified then only values in array xmins are considered when finding best xmin
* xmax: maximum value considered in calculations,values above xmax are not considered.
* seed: seed for srand function, if omitted function use srand() without parameter.

return array of simulated power law distributions with Kolmogorov smirnov test

#### bootstrap(data::AbstractArray,distribution::Type{};no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
Almost same as function above only difference is distribution parameter where should by type of distibution (con_powerlaw / dis_powerlaw). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap.

#### bootstrap(data::AbstractArray,d::UnivariateDistribution, processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
This is parallel version of bootstrap with additional parameter processes. This is number of processes which should be spawned to calculate the bootstrap. Note that this function is using [pmap](http://docs.julialang.org/en/latest/manual/parallel-computing/) which means that if there are existing workers spawned they will be used. When using this function one should consider time which is needed for data to be copied to all processes which means that it make sense to use this function only if estimate_xmin takes long time(more than few second).

#### bootstrap(data::AbstractArray,distribution::Type{},processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)  
Almost same as function above only difference is distribution parameter where should by type of distibution (con_powerlaw / dis_powerlaw). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap.

**Examples**

````jl
    f = estimate_xmin([1:100],con_powerlaw)
    bootstrap([1:100],f[1],no_of_sims = 100)
    bootstrap(collect(1:100),f[1],3,no_of_sims = 100)

    f = estimate_xmin([1:100],dis_powerlaw)
    bootstrap([1:100],f[1],no_of_sims = 100)
    bootstrap(collect(1:100),f[1],3,no_of_sims = 100)


    bootstrap(collect(1:100),con_powerlaw,3,no_of_sims = 100)
    bootstrap(collect(1:100),dis_powerlaw,no_of_sims = 100)
````

#### bootstrap_p(data::AbstractArray,d::UnivariateDistribution;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
Performs a bootstrapping hypothesis test to determine whether a power law distribution is plausible. Inspired by R [poweRlaw](http://arxiv.org/pdf/1407.3492v1.pdf) documentation.

**Parameters:**

* data: array of data which should be fit to distribution
* d: distribution, right now it support only con_powerlaw and dis_powerlaw
* no_of_sims: number of simulations which should be performed
* xmins: if not specified all unique values in data are taken as possible xmins, if specified then only values in array xmins are considered when finding best xmin
* xmax: maximum value considered in calculations,values above xmax are not considered.
* seed: seed for srand function, if omitted function use srand() without parameter.

return array of simulated power law distributions with Kolmogorov smirnov test and p_value

#### bootstrap_p(data::AbstractArray,distribution::Type{};no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
Almost same as function above only difference is distribution parameter where should by type of distibution (con_powerlaw / dis_powerlaw). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap_p.

#### bootstrap_p(data::AbstractArray,d::UnivariateDistribution, processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
This is parallel version of bootstrap_p with additional parameter processes. This is number of processes which should be spawned to calculate the bootstrap. Note that this function is using [pmap](http://docs.julialang.org/en/latest/manual/parallel-computing/) which means that if there are existing workers spawned they will be used. When using this function one should consider time which is needed for data to be copied to all processes which means that it make sense to use this function only if estimate_xmin takes long time(more than few second).

#### function bootstrap_p(data::AbstractArray,distribution::Type{},processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)  
Almost same as function above only difference is distribution parameter where should by type of distibution (con_powerlaw / dis_powerlaw). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap_p.

**Examples**

````jl
    f = estimate_xmin([1:100],dis_powerlaw)
    bootstrap_p(collect(1:100),f[1],3,no_of_sims = 100)
    bootstrap_p(collect(1:100),f[1],no_of_sims = 100)
    bootstrap_p(collect(1:100),dis_powerlaw,no_of_sims = 100)
    bootstrap_p(collect(1:100),con_powerlaw,no_of_sims = 100)
````

#### function compare_distributions(d1::ContinuousUnivariateDistribution, d2::ContinuousUnivariateDistribution, data::AbstractArray, xmin::Number = 0; sig_level = 0.05)
#### function compare_distributions(d1::DiscreteUnivariateDistribution, d2::DiscreteUnivariateDistribution, data::AbstractArray, xmin::Number = 0; sig_level = 0.05)
This function calculate Vuong test and Clarke test for non nested distributions. This is necessary since it is possible to fit power law distribution to any data set. Function was implemented according to this [Non nested model selection for spatial count regression models with application to health insurance
](https://mediatum.ub.tum.de/doc/1083601/1083601.pdf).

**Parameters:**

* d1: first distribution
* d2: second distribution
* data: data which should be used for comparing distributions
* xmin: smallest element data which was used for comparing distributions
* sig_lev: sigma level

return data type compare_distibutions with results - see type definitions in doc.

#### compare_distributions(d1::DataType, d2::DataType, data::AbstractArray, xmin::Number = 0; sig_level = 0.05)
Almost same as above functions but arguments d1 and d2 are types of distributions which are fit according to xmin and data.  

**Parameters:**

* d1: data type of first  distribution
* d2: data type of second distribution
* data: data which should be used for comparing distributions
* xmin: smallest element data which was used for comparing distributions
* sig_lev: sigma level

return data type compare_distibutions with results - see type definitions in doc.

#### function compare_distributions(d1::con_powerlaw, d2::DataType, data::AbstractArray; sig_level = 0.05)
#### function compare_distributions(d1::dis_powerlaw, d2::DataType, data::AbstractArray; sig_level = 0.05)
Almost same as above functions but argument d1 is discrete or continuous power law distribution and d2 are types of distributions which which should be compared and is fit according to xmin and data.  

**Parameters:**

* d1: first distribution (discrete or continuous power law)
* d2: data type of second distribution
* data: data which should be used for comparing distributions
* sig_lev: sigma level

return data type compare_distibutions with results - see type definitions in doc.

**Note:** distribution::Type{} means that inside curly brackets is type e.g. con_powerlaw,dis_powerlaw