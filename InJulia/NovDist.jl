"""
    NovDist

Novel distributions, copula, sampling, bootstrapping and other contents.
including:
1. Sampling
    1. LatinHyperCube(): latin hypercube sampling
2. Distributions
3. Generic methods
    1. invcdf(): inverse function of CDF, in fact, it is another name of Distributions.quantile()
		1. invstdnorm(), invnorm(), invexpo()
	2. rrank(): r-style ranking API, using ranking methosd in StatsBase module

Depends on:
1. Distributions: std lib
2. Random: std lib
3. LinearAlgebra: std lib
4. StatsBase: std lib



"""
module NovDist
    import Distributions  # std lib of distributions
	import Random  # std lib
	import LinearAlgebra  # std lib
	import StatsBase  # std lib
	# -----------
	# generic functions
	export LatinHyperCube, rrank
	# quick functions
	export invcdf, invstdnorm, invnorm, invexpo




# ----------------------- convenient methods naming
"""
	rrank( x::AbstractArray ; rev::Bool = false, ties_method::String = "ordinal" )

R-language style ranking, returns the ranks of an array's elements;
using StatsBase (std lib);
but different ties methods from R-lang;

receives:
1. x::AbstractArray, the array to get ranks
2. rev::Bool, reverse ranks? false for ordinary order, true for reverse order
3. ties_method::String, methods of ties, one of ["ordinal","dense","compete","average"]
returns:
1. Res::AbstractArray

ties methods & examples & equivalent std function:
0. example vector to rank: x = [3,4,4,99]
1. "ordinal" ==> 1,2,3,4 ==> StatsBase.ordinalrank(x) ==> "ordinal ranking"
2. "dense" ==> 1,2,2,3 ==> StatsBase.denserank(x) ==> "dense ranking"
3. "compete" ==> 1,2,2,4 ==> StatsBase.competerank(x) ==> "standard competition ranking"
4. "average" ==> 1,2.5,2.5,4 ==> StatsBase.tiedrank(x) ==> "fractional ranking"

reference:
1. http://en.wikipedia.org/wiki/Ranking#Fractional_ranking_.28.221_2.5_2.5_4.22_ranking.29
"""
function rrank( x::AbstractArray ; rev::Bool = false, ties_method::String = "ordinal" )
	if ties_method == "ordinal"
		return StatsBase.ordinalrank( x, rev = rev )::AbstractArray
	elseif ties_method == "dense"
		return StatsBase.denserank( x, rev = rev )::AbstractArray
	elseif ties_method == "compete"
		return StatsBase.competerank( x, rev = rev )::AbstractArray
	elseif ties_method == "average"
		return StatsBase.tiedrank( x, rev = rev )::AbstractArray
	else
		throw(AssertionError("invalid ties_method string"))
	end
end # end function rrank
# -----------------------







# ----------------------- Abbreviation Types
# univariate continuous distribution types
UniVarContinuousDistribution = Distributions.Distribution{Distributions.Univariate,Distributions.Continuous}

# ----------------------- Inverse CDF (Univariate distributions)
invcdf( d::Distributions.UnivariateDistribution, q::Real ) = Distributions.quantile( d, q )
# ----------------------- Quick Functions
invstdnorm( p::Real ) = invcdf( Distributions.Normal(0.0, 1.0), p )
invnorm( p::Real, μ::Real, σ::Real ) = invcdf( Distributions.Normal(μ, σ), p )
invexpo( p::Real, θ::Real ) = invcdf( Distributions.Exponential(θ), p )




# -----------------------
"""
    LatinHyperCube( D::Distributions.UnivariateDistribution, N::Int, P::Int )

Latin hypercube sampling for Univariate Distributions, e.g. Normal, Beta, Expoenential;
requires:
1. D::Distributions.UnivariateDistribution, a distribution instance
2. N::Int, a positive integer which indicates how many parts to equally divide in [0,1], the prob domain
3. P::Int, a positive integer which indicates how many samples to generate in EACH part
returns:
1. Ret::Vector{Float64}, a sampling vector with N * P sample points/elements
depends on:
1. invcdf(), of course, equivalent to Distributions.quantile()
2. Random.shuffle(), shuffles samples
timing:
1. benchmark: Intel Core i7 6700 3.40GHz, SATA3.0, no SSD, 16G DDR3
2. code: @time eval(:( D = Distributions.Normal(0.0,1.0); NovDist.LatinHyperCube(D,100,100); ) );
3. about: 0.000943 seconds (76 allocations: 160.938 KiB)
4. about: sample mu = -0.0005; sigma = 0.99803
"""
function LatinHyperCube( D::Distributions.UnivariateDistribution, N::Int, P::Int )
    # validation
    @assert( (N > 0) & (P > 0), "requries N,P > 0 but received: $(N) and $(P)" )
    # sampling in N divided Uniform(0,1)
    # NOTE: sampling N*P points in ONE Uniform(0,1) then scale them to each Uniform
    local Res::Vector{Float64} = rand( Float64, N * P )
    # loop to re-scale & inverse
    for x in 0:N-1
        for y in 1:P
			# re-scale to a new Uniform & inverse prob to a point
            Res[ x * P + y ] = Distributions.quantile( D, (Res[ x * P + y ] + x) / N )
        end
    end
	# return sampling points which are shuffled by Random.shuffle
    return Random.shuffle(Res)::Vector{Float64}
end # end function LatinHyperCube
# -----------------------
"""
	LatinHyperCube( D::Distributions.AbstractMvNormal, N::Int, P::Int )

Latin hypercube sampling for Multivariate Normal Distributions
requires:
1. D::Distributions.AbstractMvNormal, a distribution instance
2. N::Int, a positive integer which indicates how many parts to equally divide in [0,1], the prob domain, on each dimension
3. P::Int, a positive integer which indicates how many samples to generate in EACH part
returns:
1. Ret::Matrix{Float64}, a sampling vector with N * P rows and Distributions.dim(D) columns; each row is a sample point
depends on:
2. LinearAlgebra.cholesky(), Cholesky decomposition
3. rrank == StatsBase.ordinalrank(), ranks of a vector
4. LatinHyperCube( D::Distributions.UnivariateDistribution, N::Int, P::Int )
timing:
1. benchmark: Intel Core i7 6700 3.40GHz, SATA3.0, no SSD, 16G DDR3
2. code: @time eval(:( D = Distributions.MvNormal(rand(5),rand(5)); NovDist.LatinHyperCube(D,100,100); ) );
3. 17.526934 seconds (949.95 M allocations: 17.898 GiB, 6.20% gc time)
4. maximum absolute error between sampled correlation matrix and target matrix: 0.00754
reference:
1. Iman, R. L., and W. J. Conover. 1982. A Distribution-free Approach to Inducing Rank Correlation Among Input Variables. Communications in Statistics B 11:311-334
2. Zhang, Y. , & Pinder, G. . (2003). Latin hypercube lattice sample selection strategy for correlated random hydraulic conductivity fields. Water Resources Research, 39(8).
3. http://freesourcecode.net/matlabprojects/71269/latin-hypercube-sampling-in-matlab#.XLdVUuriuUk
"""
function LatinHyperCube( D::Distributions.AbstractMvNormal, N::Int, P::Int )
    # validation
	@assert( (N > 0) & (P > 0), "requries N,P > 0 but received: $(N) and $(P)" )

	# sampling in N divided Uniform(0,1)
	# NOTE: sampling N*P points in ONE Uniform(0,1) on EACH dimension
	local Ddim::Int = Distributions.dim( D.Σ )  # dim of D
	local Res::Matrix{Float64} = rand( Float64, N * P, Ddim )
	# NOTE: using Int for less memoery cost & easy indexing
	# NOTE: in practice, we do not need to declare a W0
	local Rnk1::Matrix{Int} = zeros( Int, N * P, Ddim )  # (row) ranks of each column for Res with correlation
	local EachDimStd = sqrt.( LinearAlgebra.diag( D.Σ ) )  # std of each marginal distribution

	# loop to re-scale & inverse to a MvNormal without correlation (each dim is independent)
	# NOTE: just use our LatinHyperCube() for univariate normal distribution which auto shuffles the results
	for z in 1:Ddim
		Res[:,z] = LatinHyperCube( Distributions.Normal( D.μ[z], EachDimStd[z] ), N, P )
	end

	# compute R (Ddim * Ddim size), the correlation matrix (not cov) of the Res without correlation
	# NOTE: when N*P is large, it is a good idea to use a eye(Ddim) to approximate R
	local Rmat::Matrix{Float64} = (N * P) > 10000 ?  LinearAlgebra.diagm( 0 => ones(Ddim) )  :  StatsBase.cor(Res)
	# get a new sample matrix Res1
	# NOTE: we do not specially declare Q,P but integrate them to the computing of X1
	# 		the complier is smart enough to optimize it for least memeory & time costs
	local Res1 = Res * transpose( LinearAlgebra.cholesky( Distributions.cor( D ) ).L * inv(LinearAlgebra.cholesky( Rmat ).L) )

	# record ranks of each column of Res1, then rearrange Res's columns (one by one) according to Res1's col-ranks
	for z in 1:Ddim
		# first, record X1's rank values
		Rnk1[:,z] = StatsBase.ordinalrank(Res1[:,z])
		# then, sort X0's column, the index is the rank value
		Res[:,z] = sort( Res[:,z] )
		# finally, rearrange indices according to Rnk1
		Res[:,z] = Res[Rnk1[:,z], z]
	end

	# return (we have already shuffled columns)
	return Res::Matrix{Float64}
end # end function LatinHyperCube
# -----------------------























end # end module NovDist
