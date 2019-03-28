# testing environment: Julia 1.0.1 on Intel i7-4710HQ, 12GB DDR3, Samsung SSD 860 EVO 250GB, GTX 850M

using Revise  # for debugging
push!(LOAD_PATH, "./")  # loading path

import StochInt  # import our library of: stochastic process & integral

# -------------
# DEMO: generate (discretized) daily stock returns with Brownian motion
# NOTE: usually, returns are estimated in percentages rather than digits
# TIME: 0.000044 seconds (24 allocations: 9.094 KiB)
@time StkRet = StochInt.simu_BM(
    250,   # 250 trading days
    N = 1,  # periods to discretize in one trading day, here equals one to simulate daily returns
    mu = 0.1,  # drift, estimated by percentage daily returns data
    sigma = 1.0,  # diffusion, estimated by percentage daily returns data
) .+ 2.5;  # add an initial level/returns (%) at time 0
# -------------
# DEMO: generate (discretized) daily stock prices with Geometric Brownian motion
# NOTE: in price simulation, GBM requires returns in digits but not in percentages
# TIME: 0.000018 seconds (9 allocations: 2.359 KiB)
@time StkPri = exp.( StkRet ./ 100 );
# -------------
# DEMO: generate (discretized) long-term interest rates (e.g. 10-year fed bond returns) with Ornstein-Uhlenbeck process
# NOTE: interest rates are usually estimated in percecntages rather than in digits
# TIME: 0.000046 seconds (24 allocations: 9.094 KiB)
@time IntRat = StochInt.simu_OU(
    250,  # 250 trading days
    N = 1,  # periods to discretize in one trading day
    mu = 0.1,  # drift
    sigma = 1.0,  # diffusion
) .+ 7.1;  # add benchmark interest rate level
# ------------
# # DEMO: solve the steady state distribution of a given Markov Chain (denoted by one-step transition probability matrix)
# # NOTE: element p_{i,j} in one-step matrix P means the probability from state i to state j
# # TIME:
# Pmat = [ 0.25 0.75; 0.2 0.8 ]
# @time SState = StochInt.solve_SS_MarkovChain( Pmat )
# ------------
# DEMO: solve the steady state distribution of a given Markov Chain through a fucking savage method
# NOTE: yes, just do P*P*P... until it converges
# TIME: (50 * 50 P) 0.0077 seconds (1.27 k allocations: 6.002 MiBs)
Pmat = rand(50,50); for x in 1:50; Pmat[x,:] ./= sum(Pmat[x,:]); end;  # construct sample one-step transition matrix
@time SState = StochInt.solve_SS_MarkovChain_FUCK( Pmat, atol = 1E-6, maxiter = 100 );




















#
