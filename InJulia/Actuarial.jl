"""
    Actuarial

defines basic types, methods in actuarial science (annuity, life & non-life insurances), including:
1. Types
    1. ::Annuity, standard deterministic annuity of discrete time
    2. ::ContinuousAnnuity, standard deterministic annuity of continuous time
    3. ::RandAnnuity, standard random-payment annuity of discrete time
2. Methods:
    1. AnnuityValue() & PV() & FV(), compute discounted values of a given annuity


There are some non-export members:
1. ::LifeTable, type for life tables
    1. Note: I provide data for China Life Insurance Mortality Table (2000-2003), CL(2000-2003)
    1. LT_CL1_03_NonPension::LifeTable, CL(2000-2003) for male,   for non-pension affair
    1. LT_CL2_03_NonPension::LifeTable, CL(2000-2003) for female, for non-pension affair
    1. LT_CL3_03_NonPension::LifeTable, CL(2000-2003) for male,   for pension affair
    1. LT_CL4_03_NonPension::LifeTable, CL(2000-2003) for female, for pension affair



"""
module Actuarial


# ---------------

# ---------------
"""
    Annuity

type declare of standard deterministic annuity of discrete time (pays 1 per time, start at time 0 the begin of period 1); it has fields:
1. IsFinite::Bool, if this annuity is a finite annuity, true for finite; true by default
1. IsDue::Bool  # if it is a Annuity Due (pays at the beginning of period); true for Annuity Due, false for Ordinary Annuity (pays at the end of period)
2. T::Int, the duration of this annuity, no less than 1 if IsFinite is true; when IsFinite is false, it will be neglected
3. IntRate::Float64, interest rate of this annuity, positive; discounting rate can be converted from it `` d = \\frac{i}{1+i} ``; ; for a contineous annuity, the force of interest is computed through `` \\delta = \\ln{1+i} ``

Note:
1. deferred & expired annuity can be easily computed through functions, no need to define in a type;
2. you can always construct a generic multi-inyear-payment annuity with standard annuities like constructing a portfolio!
"""
struct Annuity
    IsFinite::Bool  # type of annuity, true for finite, false for infinite
    IsDue::Bool  # if it is a Annuity Due (pays at the beginning of period); true for Annuity Due, false for Ordinary Annuity (pays at the end of period)
    T::Int  # period number, only works when IsFinite == true
    IntRate::Float64  # compound interest rate, can be converted to corresponding discounting rate `` d = \\frac{i}{1+i} ``; for a contineous annuity, the force of interest is computed through `` \\delta = \\ln{1+i} ``
    # ----
    """
        Annuity( T::Int, IntRate::Float64 ; IsFinite::Bool = true, IsDue::Bool = false )

    Constructor, where initilize a finite ordinary annuity by default;
    when IsFinite is false, T::Int won't work;
    expects a T >= 1 (when IsFinite is true), and IntRate > 0;
    expects an interest rate in digits but not in percentage, e.g. 3% -> 0.03;
    """
    function Annuity( T::Int, IntRate::Float64 ; IsFinite::Bool = true, IsDue::Bool = false )
        # validation
        IsFinite  &&  ( T >= 1 || error("Annuity type expects a T >= 1 but receives $(T)") )
        IntRate > 0.0  ||  error("Annuity type expects a positive interest rate but receives $(IntRate)")
        # new
        new( IsFinite, IsDue, T, IntRate )
    end # end constructor Annuity
end # end struct Annuity
# ---------------
"""
    AnnuityValue( a::Annuity ; AnnuityStartTime::Int = 0, DiscountTo::Int = 0 )

computes the present value of a defined annuity at a specific time/moment;

Receives:
1. a::Annuity, an Annuity instance
2. AnnuityStartTime::Int, [in Z], an integer, negative or 0 or positive, the time/moment when Annuity starts, equals 0 by default; if > 0, the annuity becomes a Deferred Annuity; if < 0, the function computes the value of this annuity at that moment
3. DiscountTo::Int, [in Z], an interger, indicating what moment to discount the given annuity to; 0 by default
Returns:
1. Res::Float64, the discounted present value of given annuity

Note:
1. this is a generic method to compute the discounted value of an annuity at any time/moment; you may easily compute a Deferred Annuity or an Expired Annuity
"""
function AnnuityValue( a::Annuity ; AnnuityStartTime::Int = 0, DiscountTo::Int = 0 )
    if a.IsFinite  # case: finite annuity
        if a.IsDue # case: finite Annuity Due
            return ( ( 1.0 - 1.0 / (1.0 + a.IntRate) ^ a.T ) / a.IntRate * (1.0 + a.IntRate)^(1.0 - AnnuityStartTime + DiscountTo) )::Float64
        else  # case: finite Ordinary Annuity
            return ( ( 1.0 - 1.0 / (1.0 + a.IntRate) ^ a.T ) / a.IntRate * (1.0 + a.IntRate)^(- AnnuityStartTime + DiscountTo) )::Float64
        end
    else  # case: infinite annuity
        if a.IsDue  # case: infinite Annuity Due
            return ( (1.0 + a.IntRate)^(1.0 - AnnuityStartTime + DiscountTo) / a.IntRate )::Float64
        else  # case: infinite Ordinary Annuity
            return ( (1.0 + a.IntRate)^(- AnnuityStartTime + DiscountTo) / a.IntRate )::Float64
        end
    end
    return -Inf::Float64  # formal return
end
# ----------------
"""
    PV( a::Annuity ; AnnuityStartTime::Int = 0 )

computes the present value of given annuity; the value is that at time 0, the beginning of period 1;

Receives:
1. a::Annuity, an Annuity instance
2. AnnuityStartTime::Int, [in Z], an integer, negative or 0 or positive, the time/moment when Annuity starts, equals 0 by default;
Returns:
1. Res::Float64, the present value of given annuity

Note:
1. you may let AnnuityStartTime > 0 to define a Deferred Annuity
"""
PV( a::Annuity ; AnnuityStartTime::Int = 0 ) = AnnuityValue( a, AnnuityStartTime = AnnuityStartTime, DiscountTo = 0 )
# ----------------
"""
    FV( a::Annuity ; AnnuityStartTime::Int = 0 )

computes the final value of given annuity; the value is that at the moment when the given annuity matures;

Receives:
1. a::Annuity, an Annuity instance
2. AnnuityStartTime::Int, [in Z], an integer, negative or 0 or positive, the time/moment when Annuity starts, equals 0 by default;
Returns:
1. Res::Float64, the present value of given annuity

Note:
1. the returned value is at the moment when the given annuity matures, not at time 0 (unless AnnuityStartTime == -1 * a.T)
2. FV() only works for finite annuity!
"""
FV( a::Annuity ; AnnuityStartTime::Int = 0 ) = a.IsFinite  ?  AnnuityValue( a, AnnuityStartTime = AnnuityStartTime, DiscountTo = AnnuityStartTime + a.T )  :  error("final value is only defined on finite annuity")
# ----------------


























end  # end module Actuarial
