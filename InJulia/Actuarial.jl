"""
    Actuarial

defines basic types, methods in actuarial science (annuity, life & non-life insurances), including:
1. Types
    1. ::CashFlow, generic type of finite cash flows, defined by pairs of (time, payment)
    2. ::InfCashFlow, generic type of infinite cash flows, equivalent to an infinite annuity
    1. ::StdAnnuity, standard deterministic annuity of discrete time
2. Methods:
    1. AnnuityValue() & PV() & FV(), compute discounted values of a given annuity


There are some non-export members:
1. ::LifeTable, type for life tables
    1. Note: I provide data for China Life Insurance Mortality Table (2000-2003), CL(2000-2003); they are experience life table for life insurance pricing
    1. LT_CL1_03_NonPension::LifeTable, CL(2000-2003) for male,   for non-pension affair
    1. LT_CL2_03_NonPension::LifeTable, CL(2000-2003) for female, for non-pension affair
    1. LT_CL3_03_NonPension::LifeTable, CL(2000-2003) for male,   for pension affair
    1. LT_CL4_03_NonPension::LifeTable, CL(2000-2003) for female, for pension affair



"""
module Actuarial
    import Base  # explicitly import Base to reload basic operators
    export cashflow, iscashflow, PV, FV, add, minus  # export generic methods
    export StdAnnuity, CashFlow  # export types

# ---------------
# """
#     LifeTable
#
# type/template for life tables, it is defined from age 0 to age Xmax,
# it has members:
# 1. Xmax::Int, the maximum age to live, equals 105 by default
# 1. q::Vector{Float64}, q_{x}, mortality in digit, from age x to age x+1, a vector with length 101 (age 0 to age 100)
# """
# struct LifeTable
#     Xmax::Int  # the maximum age to live, equals 105 by default
#     q::Vector{Float64}  # q_{x}, mortality in digit, from age x to age x+1
#     L::Vector{Float64}  # L_{x},
#
#
# end  # end struct LifeTable
# ---------------




# ---------------
CashFlow = Dict{Float64,Float64}  # define a CashFlow type, the synonym of Dict{Float64, Float64}
# ---------------
"""
    iscashflow( CF::Any )

checks if a variable is a valid CashFlow instance,
returns a Bool, true for Yes, false for No;
"""
iscashflow( CF::Any ) = begin
    isa( CF, Dict{Float64,Float64} )  ||  return false::Bool
    all( isfinite.(keys(CF)) )  ||  return false::Bool
    all( isfinite.(values(CF)) )  ||  return false::Bool
    return true::Bool
end
# ---------------
"""
    cashflow( Time::Vector{T} where T <: Real, Cash::Vector{T} where T <: Real )

external constructor of a cash flow;
in this module, we use Dict{Float64,Float64} to denote a cash flow;
requires a vector of moments (time) and corresponding cash flow values (Cash),
where duplicated time points are overwritten by the last one;
does not allow infinite time;
returns a cash flow in Dict{Float64,Float64} instance
"""
cashflow( Time::Vector{T} where T <: Real, Cash::Vector{T} where T <: Real ) = begin
    length(Time) == length(Cash)  ||  error("Time must be with the same length with Cash")
    local Res = [];
    for x in zip(Time,Cash)
        ( isfinite(x[1]) & isfinite(x[2]) )  ?  push!( Res, ( Float64(x[1]), Float64(x[2]) ) )  :  error("infinite time found in Time or Cash: $(x)")
    end
    return Dict{Float64,Float64}(Res::Vector)
end  #  end function cashflow
"""
    cashflow( OriDict::Dict{T,T} where T <: Real )

standarizes a cash flow from a dictionary
returns a cash flow in Dict{Float64,Float64} instance
"""
cashflow( OriDict::Dict{T,T} where T <: Real ) = begin
    local Res = [];
    for x in OriDict
        if isfinite(x[1]) & isfinite(x[2])
            push!(Res, ( Float64(x[2]), Float64(x[2]) ) )
        end
    end
    return Dict{Float64,Float64}(Res::Vector)
end  #  end function cashflow
# --------------
"""
    PV( CF::CashFlow, IntRate::Float64 ; DiscountTo::Union{Int,Float64} = 0.0, IsDiscrete::Bool = false )

computes the present value of a CashFlow at time CountTo::Union{Int,Float64} at a constant non-negative interest rate IntRate::Float64;

Receives:
1. CF::CashFlow, cash flow to discount
2. IntRate::Float64, a constant interest rate to discount, in digit rather than percentage
3. DiscountTo::Union{Int,Float64}, the moment to discount to
4. IsDiscrete::Bool, to discount as a discrete cash flow or a continuous cash flow, fasle by default


Note:
1. when IsDiscrete is false, money is discounted as `` PV = \\exp( - r (t-t_0) ) ``
1. when IsDiscrete is true, inner-year part is discounted at linear interpolation, e.g.
`` \\text{PV}_{t=1} = \\text{FV}_{t=1.23} \\cdot \\frac{1}{1+ 10% \\cdot 0.23 } ``
"""
function PV( CF::CashFlow, IntRate::Float64 ; DiscountTo::Union{Int,Float64} = 0.0, IsDiscrete::Bool = false )
    # validation
    iscashflow(CF)  ||  error("expects a valid CashFlow instance")
    isfinite(DiscountTo)  ||  error("expects a finite CountTo but receives $(CountTO)")
    IntRate >= 0.0  ||  error("expects a non-negative interest rate but receives: $(IntRate)")
    # declare
    local Res::Float64 = 0.0  # summed PV
    local TmpVal::Float64  # temporary value
    # compute
    if IsDiscrete  # case: discrete discounting, period by period
        for (Time,Money) in CF
            # first, discount proportional part to the floor-integer
            TmpVal = Money / ( 1.0 + IntRate * ( Time - floor(Time) ) )
            # then, discount it to floor(DiscountTo) moment
            TmpVal /= ( 1.0 + IntRate ) ^ ( floor(Time) - floor(DiscountTo) )
            # finally, discount to proportional DiscountTo moment
            Res += TmpVal * ( 1.0 + IntRate * ( DiscountTo - floor(DiscountTo) ) )
        end
    else  # case: continuous discounting, exponential
        for (Time,Money) in CF
            Res += Money * exp( IntRate * ( DiscountTo - Time ) )
        end
    end
    return Res::Float64
end  # end function PV
# --------------
"""
    add( CFa::CashFlow, CFb::CashFlow )

the addition between two CashFlow,
where cash flows at the same moment are added
"""
add( CFa::CashFlow, CFb::CashFlow ) = begin
    # no validation, auto neglect all Inf values
    local CFsum = CashFlow()
    local NewKeys::Set = union( keys(CFa), keys(CFb) )  # a union set of the two CashFlow's keys
    local tmpa::Float64, tmpb::Float64  # two temporary variables, declared here to save memories
    # compute by members
    for x in NewKeys
        tmpa = haskey(CFa, x)  ?  CFa[x]  :  0.0
        tmpb = haskey(CFb, x)  ?  CFb[x]  :  0.0
        CFsum[x] = tmpa + tmpb
    end
    return CFsum::CashFlow
end  # end function add
# --------------
"""
    minus( CFa::CashFlow, CFb::CashFlow )

the minus between two CashFlow, CFa - CFb; order matters;
"""
minus( CFa::CashFlow, CFb::CashFlow ) = begin
    # no validation, auto neglect all Inf values
    local CFminus = CashFlow()
    local NewKeys::Set = union( keys(CFa), keys(CFb) )  # a union set of the two CashFlow's keys
    local tmpa::Float64, tmpb::Float64  # two temporary variables, declared here to save memories
    # compute by members
    for x in NewKeys
        tmpa = haskey(CFa, x)  ?  CFa[x]  :  0.0
        tmpb = haskey(CFb, x)  ?  CFb[x]  :  0.0
        CFminus[x] = tmpa - tmpb
    end
    return CFminus::CashFlow
end  # end function minus
# ---------------

# ---------------
"""
    StdAnnuity

type declare of standard deterministic annuity of discrete time (pays 1 per time, start at time 0 the begin of period 1); it has fields:
1. IsFinite::Bool, if this annuity is a finite annuity, true for finite; true by default
1. IsDue::Bool  # if it is a Annuity Due (pays at the beginning of period); true for Annuity Due, false for Ordinary Annuity (pays at the end of period)
2. T::Int, the duration of this annuity, no less than 1 if IsFinite is true; when IsFinite is false, it will be neglected
3. IntRate::Float64, interest rate of this annuity, positive; discounting rate can be converted from it `` d = \\frac{i}{1+i} ``; ; for a contineous annuity, the force of interest is computed through `` \\delta = \\ln{1+i} ``

Note:
1. deferred & expired annuity can be easily computed through functions, no need to define in a type;
2. you can always construct a generic multi-inyear-payment annuity with standard annuities like constructing a portfolio!
"""
struct StdAnnuity
    IsFinite::Bool  # type of annuity, true for finite, false for infinite
    IsDue::Bool  # if it is a Annuity Due (pays at the beginning of period); true for Annuity Due, false for Ordinary Annuity (pays at the end of period)
    T::Int  # period number, only works when IsFinite == true
    IntRate::Float64  # compound interest rate, can be converted to corresponding discounting rate `` d = \\frac{i}{1+i} ``; for a contineous annuity, the force of interest is computed through `` \\delta = \\ln{1+i} ``
    # ----
    """
        StdAnnuity( T::Int, IntRate::Float64 ; IsFinite::Bool = true, IsDue::Bool = false )

    Constructor, where initilize a finite ordinary annuity by default;
    when IsFinite is false, T::Int won't work;
    expects a T >= 1 (when IsFinite is true), and IntRate > 0;
    expects an interest rate in digits but not in percentage, e.g. 3% -> 0.03;
    """
    function StdAnnuity( T::Int, IntRate::Float64 ; IsFinite::Bool = true, IsDue::Bool = false )
        # validation
        IsFinite  &&  ( T >= 1 || error("Annuity type expects a T >= 1 but receives $(T)") )
        IntRate > 0.0  ||  error("Annuity type expects a positive interest rate but receives $(IntRate)")
        # new
        new( IsFinite, IsDue, T, IntRate )
    end # end constructor StdAnnuity
end # end struct StdAnnuity
# ---------------
"""
    AnnuityValue( a::StdAnnuity ; AnnuityStartTime::Int = 0, DiscountTo::Int = 0 )

computes the present value of a defined annuity at a specific time/moment;

Receives:
1. a::StdAnnuity, an StdAnnuity instance
2. AnnuityStartTime::Int, [in Z], an integer, negative or 0 or positive, the time/moment when Annuity starts, equals 0 by default; if > 0, the annuity becomes a Deferred Annuity; if < 0, the function computes the value of this annuity at that moment
3. DiscountTo::Int, [in Z], an interger, indicating what moment to discount the given annuity to; 0 by default
Returns:
1. Res::Float64, the discounted present value of given annuity

Note:
1. this is a generic method to compute the discounted value of an annuity at any time/moment; you may easily compute a Deferred Annuity or an Expired Annuity
"""
function AnnuityValue( a::StdAnnuity ; AnnuityStartTime::Int = 0, DiscountTo::Int = 0 )
    if a.IsFinite  # case: finite annuity
        if a.IsDue # case: finite StdAnnuity Due
            return ( ( 1.0 - 1.0 / (1.0 + a.IntRate) ^ a.T ) / a.IntRate * (1.0 + a.IntRate)^(1.0 - AnnuityStartTime + DiscountTo) )::Float64
        else  # case: finite Ordinary StdAnnuity
            return ( ( 1.0 - 1.0 / (1.0 + a.IntRate) ^ a.T ) / a.IntRate * (1.0 + a.IntRate)^(- AnnuityStartTime + DiscountTo) )::Float64
        end
    else  # case: infinite annuity
        if a.IsDue  # case: infinite StdAnnuity Due
            return ( (1.0 + a.IntRate)^(1.0 - AnnuityStartTime + DiscountTo) / a.IntRate )::Float64
        else  # case: infinite Ordinary StdAnnuity
            return ( (1.0 + a.IntRate)^(- AnnuityStartTime + DiscountTo) / a.IntRate )::Float64
        end
    end
    return -Inf::Float64  # formal return
end
# ----------------
"""
    PV( a::StdAnnuity ; AnnuityStartTime::Int = 0 )

computes the present value of given annuity; the value is that at time 0, the beginning of period 1;

Receives:
1. a::StdAnnuity, an StdAnnuity instance
2. AnnuityStartTime::Int, [in Z], an integer, negative or 0 or positive, the time/moment when Annuity starts, equals 0 by default;
Returns:
1. Res::Float64, the present value of given annuity

Note:
1. you may let AnnuityStartTime > 0 to define a Deferred StdAnnuity
"""
PV( a::StdAnnuity ; AnnuityStartTime::Int = 0 ) = AnnuityValue( a, AnnuityStartTime = AnnuityStartTime, DiscountTo = 0 )
# ----------------
"""
    FV( a::StdAnnuity ; AnnuityStartTime::Int = 0 )

computes the final value of given annuity; the value is that at the moment when the given annuity matures;

Receives:
1. a::StdAnnuity, an StdAnnuity instance
2. AnnuityStartTime::Int, [in Z], an integer, negative or 0 or positive, the time/moment when Annuity starts, equals 0 by default;
Returns:
1. Res::Float64, the present value of given annuity

Note:
1. the returned value is at the moment when the given annuity matures, not at time 0 (unless AnnuityStartTime == -1 * a.T)
2. FV() only works for finite annuity!
"""
FV( a::StdAnnuity ; AnnuityStartTime::Int = 0 ) = a.IsFinite  ?  AnnuityValue( a, AnnuityStartTime = AnnuityStartTime, DiscountTo = AnnuityStartTime + a.T )  :  error("final value is only defined on finite annuity")
# ----------------


























end  # end module Actuarial
