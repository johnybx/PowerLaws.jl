using PowerLaws
using Base.Test

# tests are based moby dick word distribution
# and other data sets from http://tuvalu.santafe.edu/~aaronc/powerlaws/data.htm
# results are compared to result from R poweRlaw package
# test are done with tolerance to 1e-8
tolerance = 1e-8

moby_data = vec(readdlm(string(Pkg.dir("PowerLaws"),"/data/moby_dick.txt"),' ',Int))
cities = vec(readdlm(string(Pkg.dir("PowerLaws"),"/data/cities.txt"),' ',Int))
electrical_blackout = vec(readdlm(string(Pkg.dir("PowerLaws"),"/data/electrical_blackouts_US.txt"),' ',Int))
population = vec(readdlm(string(Pkg.dir("PowerLaws"),"/data/population.txt"),' ','\n'))

println("test: discrete powerlaw")
est = estimate_xmin(moby_data, dis_powerlaw)
@test_approx_eq_eps est[1].α 1.9527275186009831 tolerance
@test est[1].θ == 7.0
@test_approx_eq_eps est[2] 0.008252950457938835 tolerance

est = estimate_xmin(moby_data, dis_powerlaw, xmins = [2,3,4,10,20])
@test_approx_eq_eps est[1].α 1.9550379794745847 tolerance
@test est[1].θ == 10.0
@test_approx_eq_eps est[2] 0.011867106465087929 tolerance

est1 = estimate_xmin(cities, dis_powerlaw)
@test_approx_eq_eps est1[1].α 1.614392642791718 tolerance
@test est1[1].θ == 1021.0
@test_approx_eq_eps est1[2] 0.06088582981503132 tolerance

println("test: discrete powerlaw - bootstrap")
est1 = estimate_xmin(electrical_blackout, dis_powerlaw)
@test_approx_eq_eps est1[1].α 1.2201523540878356 tolerance
@test est1[1].θ == 1000.0
@test_approx_eq_eps est1[2] 0.36278306128485405 tolerance

bootstr = bootstrap(moby_data,dis_powerlaw,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,dis_powerlaw,3,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],2,no_of_sims = 15)
@test length(bootstr) == 15

println("test: continuous powerlaw")
est = estimate_xmin(moby_data, con_powerlaw)
@test_approx_eq_eps est[1].α 1.9335240463440206 tolerance
@test est[1].θ == 26.0
@test_approx_eq_eps est[2] 0.03204776085711486 tolerance

est = estimate_xmin(moby_data, con_powerlaw, xmins = [2,3,4,10,20])
@test_approx_eq_eps est[1].α 1.9514282451732778 tolerance
@test est[1].θ == 20.0
@test_approx_eq_eps est[2] 0.04416094210009813 tolerance

est = estimate_xmin(collect(1:100),con_powerlaw)
@test_approx_eq_eps est[1].α 9.178824791215408 tolerance
@test est[1].θ == 79.0
@test_approx_eq_eps est[2] 0.1824330509645657 tolerance

est1 = estimate_xmin(population,con_powerlaw)
@test_approx_eq_eps est1[1].α 2.0337071519535046 tolerance
@test est1[1].θ == 96479.0
@test_approx_eq_eps est1[2] 0.0 tolerance

println("test: continuous powerlaw - bootstrap")
bootstr = bootstrap(moby_data,con_powerlaw,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,con_powerlaw,3,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],2,no_of_sims = 15)
@test length(bootstr) == 15

println("test: compare distribution")
data = collect(1:100)
d1 = fit(dis_powerlaw,data)
d2 = fit(Poisson,data)
ll_hoods_r = logpdf(d1,data) - logpdf(d2,data)
cmpd = compare_distributions(d1,d2,data)
@test typeof(cmpd) == compare_distributions
@test cmpd.data == data
@test cmpd.log_likehoods_ratio == ll_hoods_r
@test cmpd.xmin == 1
@test cmpd.sig_level == 0.05
@test_approx_eq_eps cmpd.V_test_stat 5.7463540130003254 tolerance
@test_approx_eq_eps cmpd.V_p_val 0.9999999954405856 tolerance
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 62
@test_approx_eq_eps cmpd.C_p_val 0.020978735677851468 tolerance
@test cmpd.C_preff_distr == 1

cmpd = compare_distributions(d1,Poisson,data)
@test typeof(cmpd) == compare_distributions
@test cmpd.data == data
@test cmpd.log_likehoods_ratio == ll_hoods_r
@test cmpd.xmin == 1
@test cmpd.sig_level == 0.05
@test_approx_eq_eps cmpd.V_test_stat 5.7463540130003254 tolerance
@test_approx_eq_eps cmpd.V_p_val 0.9999999954405856 tolerance
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 62
@test_approx_eq_eps cmpd.C_p_val 0.020978735677851468 tolerance
@test cmpd.C_preff_distr == 1

cmpd = compare_distributions(dis_powerlaw,Poisson,data)
@test typeof(cmpd) == compare_distributions
@test cmpd.data == data
@test cmpd.log_likehoods_ratio == ll_hoods_r
@test cmpd.xmin == 1
@test cmpd.sig_level == 0.05
@test_approx_eq_eps cmpd.V_test_stat 5.7463540130003254 tolerance
@test_approx_eq_eps cmpd.V_p_val 0.9999999954405856 tolerance
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 62
@test_approx_eq_eps cmpd.C_p_val 0.020978735677851468 tolerance
@test cmpd.C_preff_distr == 1

cmpd = compare_distributions(dis_powerlaw,con_powerlaw,data)
@test cmpd == Union{}

cmpd = compare_distributions(d1,con_powerlaw,data)
@test cmpd == Union{}

try
    cmpd = compare_distributions(d1,con_powerlaw(),data)
catch e
    @test true == isa(e,MethodError)
end

moby_data = sort(moby_data)
d1 = estimate_xmin(moby_data,con_powerlaw)[1]
d2 = fit(Exponential,moby_data[18072:end])
cmpd = compare_distributions(d1,d2,moby_data,26)
@test cmpd.xmin == 26
@test cmpd.sig_level == 0.05
@test_approx_eq_eps cmpd.V_test_stat 8.710089666437788 tolerance
@test_approx_eq_eps cmpd.V_p_val 1.0 tolerance
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 566
@test_approx_eq_eps cmpd.C_p_val 0.0 tolerance
@test cmpd.C_preff_distr == 1

d1 = estimate_xmin(moby_data,dis_powerlaw)[1]
d2 = fit(Poisson,moby_data[15898:end])
cmpd = compare_distributions(d1,d2,moby_data,7.0)
@test cmpd.xmin == 7.0
@test cmpd.sig_level == 0.05
@test_approx_eq_eps cmpd.V_test_stat 4.448486373104684 tolerance
@test_approx_eq_eps cmpd.V_p_val 0.9999956761232807 tolerance
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 2757
@test_approx_eq_eps cmpd.C_p_val 0.0 tolerance
@test cmpd.C_preff_distr == 1
