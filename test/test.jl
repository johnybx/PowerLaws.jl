using powerlaw
using Base.Test

# tests are based moby dick word distribution
# and other data sets from http://tuvalu.santafe.edu/~aaronc/powerlaws/data.htm
# results are compared to result from R poweRlaw package
moby_data = vec(readdlm(string(Pkg.dir("powerlaw"),"/data/moby_dick.txt"),' ',Int))
cities = vec(readdlm(string(Pkg.dir("powerlaw"),"/data/cities.txt"),' ',Int))
electrical_blackout = vec(readdlm(string(Pkg.dir("powerlaw"),"/data/electrical_blackouts_US.txt"),' ',Int))
population = vec(readdlm(string(Pkg.dir("powerlaw"),"/data/population.txt"),' ','\n'))

println("test: discrete powerlaw")
est = estimate_xmin(moby_data, dis_powerlaw)
@test est[1].α == 1.9527275186009831
@test est[1].θ == 7.0
@test est[2] == 0.008252950457938835

est = estimate_xmin(moby_data, dis_powerlaw, xmins = [2,3,4,10,20])
@test est[1].α == 1.9550379794745847
@test est[1].θ == 10.0
@test est[2] == 0.011867106465087929

est1 = estimate_xmin(cities, dis_powerlaw)
@test est1[1].α == 1.614392642791718
@test est1[1].θ == 1021.0
@test est1[2] == 0.06088582981503132

println("test: discrete powerlaw - bootstrap")
est1 = estimate_xmin(electrical_blackout, dis_powerlaw)
@test est1[1].α == 1.2201523540878356
@test est1[1].θ == 1000.0
@test est1[2] == 0.36278306128485405

bootstr = bootstrap(moby_data,dis_powerlaw,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,dis_powerlaw,3,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],2,no_of_sims = 15)
@test length(bootstr) == 15

println("test continuous powerlaw")
est = estimate_xmin(moby_data, con_powerlaw)
@test est[1].α == 1.9335240463440206
@test est[1].θ == 26.0
@test est[2] == 0.03204776085711486

est = estimate_xmin(moby_data, con_powerlaw, xmins = [2,3,4,10,20])
@test est[1].α == 1.951428237468947
@test est[1].θ == 20.0
@test est[2] == 0.04416094210009813

est = estimate_xmin(collect(1:100),con_powerlaw)
@test est[1].α == 9.178824791215408
@test est[1].θ == 79.0
@test est[2] == 0.1824330509645657

est1 = estimate_xmin(population,con_powerlaw)
@test est1[1].α == 2.0337071519535046
@test est1[1].θ == 96479.0
@test est1[2] == 0.0

println("test continuous powerlaw - bootstrap")
bootstr = bootstrap(moby_data,con_powerlaw,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,con_powerlaw,3,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],2,no_of_sims = 15)
@test length(bootstr) == 15
