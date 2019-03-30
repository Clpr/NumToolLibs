"""
    StochInt

StochInt is a collection of functions about stochastic integral. It consists of:
1. simu_BM(): generate samples of Brownian motion (discretized) which is a drift-diffusion process
2. simu_OU(): generate samples of Ornstein-Uhlenbeck process (discretized) which is a mean-reverting process
3. solve_SS_MarkovChain(): solve steady state (stationary distribution) for a given Markov chain (of course, discrete)
4. solve_SS_MarkovChain_FUCK(): a fucking savage algorithm to compute steady state, through computing P^n for enough n degrees until a tolerance is met
5. simu_MarkovChain(): generate samples of a Markov Chain, requires initial state (integer from 1,...,k; where k is the number of all states )

And there are some non-export underlying functions:
1. randDiscrete(): generate a random integer from a given discrete distribution, e.g. [0.1, 0.2, 0.4, 0.3] --> 2 from {1,2,3,4}


There are some notes:
1. I do not provide functions for Geometric Brownian Motion (GBM), because it can always be converted from a Brownian motion by exponentialization;



References:
1. "Stokey, N. L. . (2008). The economics of inaction: stochastic control models with fixed costs."
2. "Press, William H., Flannery, B. P., Teukolsky, S. A., & Vetterling, W. T. (1988). Numerical recipes in C: the art of scientific computing. Numerical Recipes in C The Art of Scientific Computing."


If any question or issue, please email to:
{tianhao.zhao@emory.edu}, Tianhao Zhao (GitHub: Clpr)
"""
module StochInt


# ----------
"""
    randDiscrete( pVec::Vector{Float64} )

generates a random integer based on a given discrete distribution.
e.g. input a pVec = [0.3, 0.4, 0.3], select one element from {1,2,3},
where {1,2,3} correspond to probabilities [0.3,0.4,0.3] in order.

Receives:
1. pVec::Vector{Float64}, discrete distribution function
Returns:
1. Res::Int, a random integer

Timing:
about 0.0015 seconds (10.00 k allocations: 1.068 MiB) for F = [0.1,0.2,0.3,0.4] per 10000 times;
the testing is in main namespace, which is much slower than called by other functions in practice.


Mathematics:
for a discrete distribution like F = {1,...,k} -> [p1,...,pk] where sum(F) = 1.0,
firstly generate a uniform random number R from [0,1] uniform distribution,
then try to locate R between two probabilities pi & pj;
"""
randDiscrete( pVec::Vector{Float64} ) = begin
    # validation
    @assert( isapprox(1.0, sum(pVec)) , "requires a discrete distribution whose added probabilities is 1, but received $(sum(pVec))" )
    # generate
    local kState = length(pVec)
    local pVecSum = insert!( cumsum(pVec), 1, -1E-8 )  # convert to CDF & add a first minor negative number (to avoid the branch for exact 0)
    # to avoid the branch for exact 1
    pVecSum[end] += 1E-8
    local RndUnif = rand()  # generate a random number from [0,1] uniform distribution
    local Res = 1::Int # result

    while Res < kState  # totally compare T times (including in [0.0, p1])
        if pVecSum[Res] <= RndUnif < pVecSum[Res + 1]
            break
        end
        Res += 1
    end
    return Res::Int
end # end function randDiscrete
# ----------
"""
    simu_BM( T::Int ; N::Int = 1, mu::Float64 = 0.0, sigma::Float64 = 1.0 )

approximates a Brownian motion `` dX_t = \\mu dt + \\sigma dW_t `` with parameter ``( \\mu, \\sigma^2 )``,
where ``X_t`` is a Brownian motion with **constant** drift ``\\mu`` and **constant** diffusion ``\\sigma``, ``W_t`` is a standard Gaussian distribution;

Receives:
1. T::Int, how many periods to simulate (e.g. years)
2. N::Int, how many points during one period, by default 1
3. mu::Float64, drift term, a constant, by default 0.0
4. sigma::Float64, diffusion term, a constant, by default 1.0
Returns:
1. Vector{Float64} = Array{Float64,1}, a simulated path, with length T * N, starting from X0 = 0

According to (Stokey, 2008) Ch 2.4, I use ``X_{t+\\Delta t} = X_{t} + \\Delta X`` approximate the defined Brownian motion,
where ``\\Delta X = \\begin{cases}  h, P(\\Delta X = h) = p  \\\\  -h, P(\\Delta X = -h) = 1 - p     \\end{cases}`` is a random variable.
We have ``h = \\begin{cases} \\frac{1}{2}( 1 + \\frac{\\mu\\sqrt{\\Delta t}}{\\sqrt{ \\sigma^2 + \\mu^2 \\Delta t }} )   \\\\  \\approx  \\frac{1}{2}( 1 + \\frac{\\mu}{\\sigma}\\sqrt{\\Delta t} )  ,  \\Delta t \\ll \\frac{\\sigma^2}{\\mu^2}  \\end{cases}``,
and `` h = \\sigma \\sqrt{\\Delta t} ``;
"""
function simu_BM( T::Int ; N::Int = 1, mu::Float64 = 0.0, sigma::Float64 = 1.0 )
    # decalre
    local DeltaT = 1 / N  # step length
    local Res = [0.0]  # result, use push!() then to acquire higher performance
    local RndShock = rand(Float64, T * N - 1 )  # generating random numbers first; NOTE: minus 1 because it is increments
    # NOTE: I do not use conventional approximation ``pThreshold = \\frac{1}{2}()`` to allow user-friendly simulations of daily stock returns which has DeltaT = 1 (usually larger than ``\\sigma^2/\\mu^2``)
    local pThreshold = 0.5 * ( 1.0 + mu * sqrt(DeltaT) / sqrt( sigma^2 + mu^2 * DeltaT) )  # probability threshold to RISE hMove::Float64, for each discretized delta-t period
    local hRise = sigma * sqrt(DeltaT)  # movement size, if pThreshold happens, rise hRise; otherwise, fall hRise
    # generate Delta X
    for t in 1:(T * N - 1)
        if RndShock[t] <= pThreshold
            push!( Res, Res[t] + hRise )  # rise case
        else
            push!( Res, Res[t] - hRise )  # fall case
        end
    end
    # return
    return Res::Vector{Float64}
end  # end function simu_BM
# -----------
"""
    simu_OU( T::Int ; N::Int = 1, mu::Float64 = 0.0, sigma::Float64 = 1.0 )

generates simulated series following a Ornstein-Uhlenbeck process with parameter ``(\\mu, \\sigma^2)``,
where ``\\mu`` is drift term, ``\\sigma`` is diffusion term;
According to (Stokey, 2008), we use algorithms similar to that in simu_BM()::Function to approximate
`` dX_t = -\\mu dt + \\sigma dW_t ``, where ``W_t`` is a Wiener process.

Receives:
1. T::Int, how many periods to simulate (e.g. years)
2. N::Int, how many points during one period, by default 1
3. mu::Float64, drift term, a constant, by default 0.0
4. sigma::Float64, diffusion term, a constant, by default 1.0
Returns:
1. Vector{Float64} = Array{Float64,1}, a simulated path, with length T * N, starting from X0 = 0

Mathematics:
1. ``h = \\sigma \\sqrt{\\Delta t}``
2. `` p(X_t) = 0.5 - \\frac{\\alpha X_t}{ 2\\sigma } \\sqrt{\\Delta t} ``
3. `` \\Delta X_t = \\begin{cases} h, \\text{Prob}=p(X_t)  \\\\ -h, \\text{Prob} = 1 - p(X_t)  \\end{cases} ``
3. then, `` X_{t+\\Delta t} = X_t + \\Delta X_t ``

"""
function simu_OU( T::Int ; N::Int = 1, mu::Float64 = 0.0, sigma::Float64 = 1.0 )
    # decalre
    local DeltaT = 1 / N  # step length
    local Res = [0.0]  # result, use push!() then to acquire higher performance
    local RndShock = rand(Float64, T * N - 1 )  # generating random numbers first; NOTE: minus 1 because it is increments
    # NOTE: different from BM & GBM, the pThreshold::Float64 for a OU process may depend on X_t, which means we have to dynamically compute it
    local hRise = sigma * sqrt(DeltaT)  # movement size, if pThreshold happens, rise hRise; otherwise, fall hRise
    # generate Delta X
    for t in 1:( T * N - 1 )
        if RndShock[t] <= (  0.5 - mu * Res[t] * sqrt(DeltaT) / 2 / sigma  )  # this is pThreshold
            push!( Res, Res[t] + hRise )
        else
            push!( Res, Res[t] - hRise )
        end
    end
    return Res::Vector{Float64}
end  # end function simu_OU
# -----------------
# NOTE: need debug (math is wrong)
# """
#     solve_SS_MarkovChain( Pmat::Array{Float64,2} )
#
# solves the steady state (stationary distribution) for a given Markov chain,
# where the Markov chain is denoted with a one-step matrix ``P``.
# In matrix P, element x{row,coloumn} indicates the one-step transition probability
# from state row to state column;
#
# Receives:
# 1. Pmat::Array{T,2} where T, one-step transition probability matrix
# Returns:
# 1. Res::Vector{Float64}, vector of steady state (probabilities), states are sorted along with the order of matrix P's rows
#
# Notes:
# 1. we just solve the equations `` \\vec{\\pi} \\cdot \\bm{P} = \\vec{\\pi} ``
# 2. following conventions, we add a new all-one column to matrix P to constrain the steady state: `` \\sum^k \\pi_k = 1 ``
# 2. only least validations (e.g. inversibility of matrix P) to acquire better performance
# 3. without loss of generality, we only provide a Float64 version function
#
# Mathematics (how to solve `` \\vec{\\pi} \\cdot \\bm{P} = \\vec{\\pi} ``):
# 1. write element-level equation for a specific state j in total k states: `` \\pi_1 p_{1,j} + \\dots + \\pi_jp_{j,j} + \\dots + \\pi_{k} p_{k,j} = \\pi_j ``
# 2. move ``\\pi_j`` to the left side: `` \\pi_1 p_{1,j} + \\dots + \\pi_j(p_{j,j}-1) + \\dots + \\pi_{k} p_{k,j} = 0  ``
# 3. substitute the constraint of ``\\pi_j``: ``\\sum^k_j \\pi_j = 1``
# 4. get: `` \\pi_1 p_{1,j} + \\dots + (1-\\sum^{i\\neq j} \\pi_i)(p_{j,j}-1) + \\dots + \\pi_{k} p_{k,j} = 0   ``
# 5. simplify and get: `` \\pi_1 (p_{1,j} - p_{j,j} + 1) + \\dots + \\pi_{j} 0 + \\dots + \\pi_k(p_{k,j} - p_{j,j} + 1) = 1 - p_{j,j}  ``
# 6. write the equations in a linear algebra form: `` \\pi_{1\\times k} A_{k \\times k} = b_{1 \\times k} ``
# 7. to compute matrix A, firtst compute a matrix A': `` A'_{k\\times k} = P_{k\\times k} + 1_{k\\times k} - 1_{k\\times 1}\\cdot \\text{diag}(P_{k\\times k})_{1\\times k} ``
# 8. then matrix A is: `` A_{k\\times k} = A'_{k\\times k} - \\text{diag}(A')_{k\\times k} ``
# 9. matrix b is: `` b_{1\\times k} = 1_{1\\times k} - \\text{diag}(P_{k\\times k})_{1\\times k} ``
# 10. please note, index j traversals on rows, and index k traversals on columns
# """
# function solve_SS_MarkovChain( Pmat::Array{Float64,2} )
#     # declare
#     local kState = size(Pmat)  # number of states
#     # assertion
#     @assert( kState[1] == kState[2], " requires an n*n Array{Float64,2} matrix but received a $(kState) $(typeof(Pmat)) " )
#     # construct matrix A & matrix b
#     local Amat = copy(Pmat); local bmat = zeros(Float64, 1, kState[1])
#     for idxj in 1:kState[1]  # loop on rows/states
#         # construct matrix b
#         bmat[1,idxj] = 1.0 - Amat[idxj,idxj]
#         # construct matrix A'
#         for idxk in 1:kState[2]  # loop on columns (though A is a n*n matrix, we use [1],[2] to distinguish rows & columns)
#             Amat[idxj,idxk] += 1.0 - Amat[idxj,idxj]
#         end
#         # then set diag(A) as 0 (equals to A = A' - diag(A'))
#         Amat[idxj,idxj] = 0.0
#     end
#     # solve the linear equation system
#     local Res::Vector{Float64} = ( bmat * inv(Amat) )[1,:]  # squeeze to 1-D vector
#
#     return Res::Vector{Float64}
# end # end function solve_SS_MarkovChain
# --------------------
"""
    solve_SS_MarkovChain_FUCK( Pmat::Array{Float64,2} ; atol::Float64 = 1E-6, maxiter::Int = 100 )

solves the steady state of a given Markov chain with one-step transition probability matrix Pmat::Array{Float64,2},
where element p_{i,j} denotes the probability from state i to state j.
uses a fucking savage algorithm: computing P*P*P.... until converges.

Receives:
1. Pmat::Array{Float64,2}, one-step transition probability matrix
2. atol::Float64, tolerance to converge
3. maxiter::Int, maximum iteration times to stop
Retruns:
1. Res::Vector{Float64}, result vector (steady state distribution)

(Actually, fuckingly naive but very simple & high performance!)
"""
function solve_SS_MarkovChain_FUCK( Pmat::Array{Float64,2} ; atol::Float64 = 1E-6, maxiter::Int = 100 )
    # assertion (one-step transition matrix definition)
    local kState = size(Pmat)  # number of states
    local Pmat2 = copy(Pmat)  # temp variable
    @assert( kState[1] == kState[2], " requires an n*n Array{Float64,2} matrix but received a $(kState) $(typeof(Pmat)) " )
    @assert( all( abs.(sum(Pmat, dims=2) .- 1.0) .< 1E-6 ), "not a one-step transition matrix, requires the sum of every row is one!")
    # decalre
    local Res = zeros(Float64,kState[1]);  # result, also temp variable
    local flag_exit = true;  # convergence flag, false to stop loops
    local counter = 1;  # counter
    # fuckingly compute
    while flag_exit
        Res[:] .= Pmat2[1,:]  # current value
        Pmat2[:,:] *= Pmat2 * Pmat  # update
        if counter >= maxiter | isapprox(Pmat2[1,:], Res, atol = atol)  # if converged
            flag_exit = false
        end
        counter += 1
    end
    return Res::Vector{Float64}
end  # end function solve_SS_MarkovChain_FUCK
# ----------------
"""
    simu_MarkovChain( Pmat::Matrix{Float64}, InitState::Int, T::Int )

generates a random Markov Chain when given one-step transition matrix, an initial state and periods of time to simulate.

Receives:
1. Pmat::Matrix{Float64}, one-step transition probability matrix, where P{i,j} means the probability from state i to state j
2. InitState::Int, initial state to start, from 1 to size(Pmat)[1] (row/column number of Pmat)
3. T::Int, time periods to simulate, a finite positive integer
Returns:
1. Res::Vector{Int}, a simulated Markov Chain with length T
Dependency [no external]:
1. StochInt.randDiscrete(), generate a random integer from a given discrete distribution

"""
function simu_MarkovChain( Pmat::Matrix{Float64}, InitState::Int, T::Int )
    # declare
    local kState = size(Pmat)  # size of one-step transition matrix
    local Res = [InitState]  # result, a Vector{Int}
    # validation (one-step transition matrix definition)
    @assert( kState[1] == kState[2], " requires an n*n Array{Float64,2} matrix but received a $(kState) $(typeof(Pmat)) " )
    @assert( all( abs.(sum(Pmat, dims=2) .- 1.0) .< 1E-6 ), "not a one-step transition matrix, requires the sum of every row is one!")
    # generate
    for t in 1:T-1 # do not use 2:T to avoid bounding error
        # select a discrete distribution from the transition matrix according to current state
        # then randomly generate a new state
        push!(Res, randDiscrete( Pmat[Res[end],:] ) )  # [:] auto squeezes the row-vector to a one-dimension vector
    end
    return Res::Vector{Int}
end # end function simu_MarkovChain
# -------------
















end  # end module: StochInt
