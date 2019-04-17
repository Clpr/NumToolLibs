"""
    NovDist

Novel distributions, copula, sampling, bootstrapping and other contents.
including:
1. Sampling
    1. LatinHyperCube(): latin hypercube sampling
2. Distributions
3. Generic methods
    1. invcdf(): inverse function of CDF, in fact, it is another name of Distributions.quantile()


Depends on:
1. Distributions: std lib



"""
module NovDist
    import Distributions  # std lib of distributions
	# -----------
	# generic functions
	export invcdf, LatinHyperCube
	# quick functions
	export invstdnorm, invnorm, invexpo



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

Latin hypercube sampling for Univariate Continuous Distributions, e.g. Normal, Beta, Expoenential;
requires:
1. D::Distributions.UnivariateDistribution, a distribution instance
2. N::Int, a positive integer which indicates how many parts to equally divide in [0,1], the prob domain
3. P::Int, a positive integer which indicates how many samples to generate in EACH part
returns:
1. Ret::Vector{Float64}, a sampling vector with N * P sample points/elements
depends on:
1. invcdf(), of course, equivalent to Distributions.quantile()
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
	# return sampling points
    return Res::Vector{Float64}
end # end function LatinHyperCube
# -----------------------































end # end module NovDist
