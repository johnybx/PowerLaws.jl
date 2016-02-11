using powerlaw
using Base.Test

# tests are based moby dick word distribution
# results are compared to result from R poweRlaw package
moby_data = vec(readdlm(string(Pkg.dir("powerlaw"),"/data/moby_dick.txt"),' ','\n'))
est = estimate_xmin(moby_data,dis_powerlaw)
@test est[1].α == 1.9527275186009831 
@test est[1].θ == 7.0
@test est[2] == 0.008252950457938835
