TimeCost = []
D = Distributions.Normal(0.0, 1.0)
for x in 1:10000
  # sampling, and get the time cost (gc time included), then store
  push!(TimeCost, @timed( NovDist.LatinHyperCube(D, 100,100) )[2] )
end
# averaging
StatsBase.mean(TimeCost)



PyPlot.plt[:hist]( tmp, 40, density = true );
PyPlot.grid(true); PyPlot.ylabel("Density"); PyPlot.xlabel("x");
PyPlot.title("PDF of 10000 samples from N(0,1)");
PyPlot.tight_layout()



TimeCost = []
Sigma = [0.5 -0.25; -0.25 1.5];
Mu = [ 0.0, 1.0 ];
D = Distributions.MvNormal(Mu, Sigma)
for x in 1:10000
  push!(TimeCost, @timed( NovDist.LatinHyperCube(D, 100,100) )[2] )
end
StatsBase.mean(TimeCost)
StatsBase.std(TimeCost)










#
