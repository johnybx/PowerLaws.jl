Documentation
========

### Type definition

This definitions are done according to [Distributions](http://distributionsjl.readthedocs.org/en/latest/) package template. Continuous powerlaw support all function like cdf,ccdf,pdf,median,var etc. Discrete powerlaw do not have implemented moments - median,variance etc.

#### con_powerlaw(alfa,xmin)
Continuios power law with probability density function ![equation](http://www.sciweavers.org/tex2img.php?eq=p%28x%29%20%3D\begin{cases}%20\frac{%28\alpha%20%20-%201%29}{\theta}%20*%20%28\frac{x}{\theta}%29^{-\alpha}%20%20%20%26%20x%20%20\geq%20%20\theta%20\\0%20%26%20x%20%20%3C%20\theta%20\end{cases}%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0%22%20align=%22center%22%20border=%220%22%20alt=%22p%28x%29%20=\begin{cases}%20\frac{%28\alpha%20%20-%201%29}{\theta}%20*%20%28\frac{x}{\theta}%29^{-\alpha}%20%20%20&%20x%20%20\geq%20%20\theta%20\\0%20&%20x%20%20%3C%20\theta%20\end{cases}%20%22%20width=%22244%22%20height=%2250%22)

#### dis_powerlaw(alfa,xmin)
Discrete power law with probability mass function ![equation](http://www.sciweavers.org/tex2img.php?eq=p%28x%29%20%3D\begin{cases}%20\frac{x^{-\alpha}}{zeta%28%20\alpha%20%2C%20\theta%20%29}%20%20%20%26%20x%20%20\geq%20%20\theta%20\\0%20%26%20x%20%20%3C%20\theta%20\end{cases}%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0%22%20align=%22center%22%20border=%220%22%20alt=%22p%28x%29%20=\begin{cases}%20\frac{x^{-\alpha}}{zeta%28%20\alpha%20,%20\theta%20%29}%20%20%20&%20x%20%20\geq%20%20\theta%20\\0%20&%20x%20%20%3C%20\theta%20\end{cases}%20%22%20width=%22201%22%20height=%2254%22)

### Functions

#### estimate_xmin(data::AbstractArray,discrete::Bool = false,xmins::AbstractArray = [],xmax::Int64 = int(1e5))

Estimates best xmin and alfa for data set with respect to the [Kolmogorov smirnov test](https://www.encyclopediaofmath.org/index.php/Kolmogorov-Smirnov_test)

**Parameters:**

* data: array of data which should be fit to distribution
* discrete: if true data are considered discrete, continuous otherwise
* xmins: if not specified all unique values in data are taken as possible xmins, if specified then only values in array xmins are considered when finding best xmin
* xmax: maximum value considered in calculations,values above xmax are not considered in for example calculating Kolmogorov smirnov test.

return powerlaw distribution (continuos or discrete) with best xmin and alfa and Kolmogorov smirnov test which belongs to the distribution.
    
#### Kolmogorov_smirnov_test{T<:UnivariateDistribution}(dat::AbstractArray,d::T,xmax::Int64 = int(1e5))

Calculate [Kolmogorov smirnov test](https://www.encyclopediaofmath.org/index.php/Kolmogorov-Smirnov_test) on given data and distribution.

**Parameters:**

* data: array of data which should be fit to distribution
* d: distribution, right now it support only con_powerlaw and dis_powerlaw
* xmax: maximum value considered in calculations,values above xmax are not considered.

return KS_statistic

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

#### function bootstrap(data::AbstractArray,discrete::Bool;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
Almost same as function above only difference is that instead of distribution parameter there is discreate parameter which indicates whether data are discrete or continuous (true for discrete and false for continuous). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap.

#### bootstrap(data::AbstractArray,d::UnivariateDistribution, processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
This is parallel version of bootstrap with additional parameter processes. This is number of processes which should be spawned to calculate the bootstrap. Note that this function is using [pmap](http://docs.julialang.org/en/latest/manual/parallel-computing/) which means that if there are existing workers spawned they will be used. When using this function one should consider time which is needed for data to be copied to all processes which means that it make sense to use this function only if estimate_xmin takes long time(more than few second).

#### function bootstrap(data::AbstractArray,discrete::Bool,processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)  
Almost same as function above only difference is that instead of distribution parameter there is discreate parameter which indicates whether data are discrete or continuous (true for discrete and false for continuous). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap.


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

#### function bootstrap_p(data::AbstractArray,discrete::Bool;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
Almost same as function above only difference is that instead of distribution parameter there is discreate parameter which indicates whether data are discrete or continuous (true for discrete and false for continuous). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap_p.

#### bootstrap_p(data::AbstractArray,d::UnivariateDistribution, processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)
This is parallel version of bootstrap_p with additional parameter processes. This is number of processes which should be spawned to calculate the bootstrap. Note that this function is using [pmap](http://docs.julialang.org/en/latest/manual/parallel-computing/) which means that if there are existing workers spawned they will be used. When using this function one should consider time which is needed for data to be copied to all processes which means that it make sense to use this function only if estimate_xmin takes long time(more than few second).

#### function bootstrap_p(data::AbstractArray,discrete::Bool,processes::Int64;no_of_sims::Int64 = 10,xmins::AbstractArray = [],xmax::Int64 = round(Int,1e5),seed::Int64 = 0)  
Almost same as function above only difference is that instead of distribution parameter there is discreate parameter which indicates whether data are discrete or continuous (true for discrete and false for continuous). This function calls estimate_xmin on data with parameter distrete and than perform bootstrap_p.

