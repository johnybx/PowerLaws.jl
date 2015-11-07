Documentation
========

### Type definition

This definitions are done according to [Distributions](http://distributionsjl.readthedocs.org/en/latest/) package template. Continuous powerlaw support all function like cdf,ccdf,pdf,median,var etc. Discrete powerlaw do not have implemented moments - median,variance etc.

#### con_powerlaw(alfa,xmin)
Continuios power law with probability density function [equation](http://www.sciweavers.org/tex2img.php?eq=p%28x%29%20%3D\begin{cases}%20\frac{%28\alpha%20%20-%201%29}{\theta}%20*%20%28\frac{x}{\theta}%29^{-\alpha}%20%20%20%26%20x%20%20\geq%20%20\theta%20\\0%20%26%20x%20%20%3C%20\theta%20\end{cases}%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0%22%20align=%22center%22%20border=%220%22%20alt=%22p%28x%29%20=\begin{cases}%20\frac{%28\alpha%20%20-%201%29}{\theta}%20*%20%28\frac{x}{\theta}%29^{-\alpha}%20%20%20&%20x%20%20\geq%20%20\theta%20\\0%20&%20x%20%20%3C%20\theta%20\end{cases}%20%22%20width=%22244%22%20height=%2250%22)

#### dis_powerlaw(alfa,xmin)
Discrete power law with probability mass function ![equation](http://www.sciweavers.org/tex2img.php?eq=p%28x%29%20%3D%5Cbegin%7Bcases%7D%20%5Cfrac%7Bx%5E%7B-%5Calpha%7D%7D%7Bzeta%28%20%5Calpha%20%2C%20%5Ctheta%20%29%7D%20%20%20%26%20x%20%20%5Cgeq%20%20%5Ctheta%20%5C%5C0%20%26%20x%20%20%3C%20%5Ctheta%20%5Cend%7Bcases%7D%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0%22%20align=%22center%22%20border=%220%22%20alt=%22p(x)%20=\begin{cases}%20\frac{x^{-\alpha}}{zeta(%20\alpha%20,%20\theta%20)}%20%20%20&%20x%20%20\geq%20%20\theta%20\\0%20&%20x%20%20%3C%20\theta%20\end{cases}%20%22%20width=%22201%22%20height=%2254%22)

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

#### function bootstrap(data::AbstractArray,distribution::Type{};no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
Almost same as function above only difference is that instead of distribution parameter there is discreate parameter which indicates whether data are discrete or continuous (true for discrete and false for continuous). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap.

#### bootstrap(data::AbstractArray,d::UnivariateDistribution, processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
This is parallel version of bootstrap with additional parameter processes. This is number of processes which should be spawned to calculate the bootstrap. Note that this function is using [pmap](http://docs.julialang.org/en/latest/manual/parallel-computing/) which means that if there are existing workers spawned they will be used. When using this function one should consider time which is needed for data to be copied to all processes which means that it make sense to use this function only if estimate_xmin takes long time(more than few second).

#### function bootstrap(data::AbstractArray,distribution::Type{},processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)  
Almost same as function above only difference is that instead of distribution parameter there is discreate parameter which indicates whether data are discrete or continuous (true for discrete and false for continuous). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap.

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

#### function bootstrap_p(data::AbstractArray,distribution::Type{};no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
Almost same as function above only difference is that instead of distribution parameter there is discreate parameter which indicates whether data are discrete or continuous (true for discrete and false for continuous). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap_p.

#### bootstrap_p(data::AbstractArray,d::UnivariateDistribution, processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
This is parallel version of bootstrap_p with additional parameter processes. This is number of processes which should be spawned to calculate the bootstrap. Note that this function is using [pmap](http://docs.julialang.org/en/latest/manual/parallel-computing/) which means that if there are existing workers spawned they will be used. When using this function one should consider time which is needed for data to be copied to all processes which means that it make sense to use this function only if estimate_xmin takes long time(more than few second).

#### function bootstrap_p(data::AbstractArray,distribution::Type{},processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)  
Almost same as function above only difference is that instead of distribution parameter there is discreate parameter which indicates whether data are discrete or continuous (true for discrete and false for continuous). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap_p.

**Examples**

````jl
    f = estimate_xmin([1:100],dis_powerlaw)
    bootstrap_p(collect(1:100),f[1],3,no_of_sims = 100)
    bootstrap_p(collect(1:100),f[1],no_of_sims = 100)
    bootstrap_p(collect(1:100),dis_powerlaw,no_of_sims = 100)
    bootstrap_p(collect(1:100),con_powerlaw,no_of_sims = 100)
````


**Note:** distribution::Type{} means that inside curly brackets is type e.g. con_powerlaw,dis_powerlaw