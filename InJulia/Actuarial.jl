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
    import Base  # explicitly import Base to reload basic operators (such as +, -)
    export cashflow, iscashflow, PV, FV, cfcompress  # export generic methods
    export StdAnnuity, Payment, CashFlow  # export types

# ==================================
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




# ==================================
"""
    Payment( Time::Real, Payment::Real )

a synonym of Pair{Float64, Float64}, a "Time => Payment" relationship;
the constructor can automatically convert input args to Float64;
allows finite/infinite Time/Payment; please refer to the doc of CashFlow for more information;
"""
Payment = Pair{Float64, Float64}
# --------------- Atomic Operations
# NOTE: overloads basic Arithmetical operators for Payment
# NOTE: two payments cannot be added, so, no def for + & -
Base.:-( pay::Payment ) = Payment( pay.first, - pay.second )  # but it can be negative!
Base.:*( pay::Payment, k::Real ) = Payment( pay.first, Float64(k * pay.second) )
Base.:*( k::Real, pay::Payment ) = Payment( pay.first, Float64(k * pay.second) )  # exchange-able principle does not always stand, needs to define it manually
Base.:/( pay::Payment, k::Real ) = Payment( pay.first, Float64(pay.second / k) )
Base.:/( k::Real, pay::Payment ) = Payment( pay.first, Float64(k / pay.second) )
"""
    PV( CF::CashFlow, IntRate::Float64 ; DiscountTo::Union{Int,Float64} = 0.0, IsDiscrete::Bool = false )

computes given payment's present value at a finite moment DiscountTo::Float64;

Receives:
1. CF::Payment, a given payment
2. IntRate::Float64, a constant interest rate to discount at, >= 0
3. DiscountTo::Real, the moment to discount to, 0.0 by default
4. IsDiscrete::Bool, how to deal with the non-integer (digit) part of period, true by default
Returns:
1. Res::Float64, PV of the given payment

Mathematics:
1. if IsDiscrete is false, uses exponential discounting: `` {PV} = \\exp(-rt) ``, where r is IntRate
2. if IsDiscrete is true, uses linear interpolation to discount inner-period (non-integer, digit) part: `` {PV}_{t=1} = {PV}_{t=1.25} / ( 1 + r * 0.25 ) ``, where r is IntRate
3. for infinite time/payment, we have principles (we use (time,payment) to denote a Payment):
    1. we use R to denote finite real number
    1. non-zero R = finite DiscountTo + one of: (R,R)
    2. zero = finite DiscountTo + one of: (Inf,R)
    3. Inf = finite DiscountTo + one of: (-Inf,Inf), (R,Inf)
    5. -Inf = finite DiscountTo + one of: (R,-Inf),(-Inf,-Inf)
    6. NaN = finite DiscountTo + one of: (Inf,Inf),(Inf,-Inf),(-Inf,R)
"""
PV( CF::Payment, IntRate::Float64 ; DiscountTo::Real = 0.0, IsDiscrete::Bool = true ) = begin
    # valid & declare
    isfinite(DiscountTo)  ||  error("expects a finite DiscountTo but received: $(DiscountTo)")
    IntRate >= 0.0  ||  error("expects a non-negative interest rate but received: $(IntRate)")
    # compute
    if isfinite(CF.first) & isfinite(CF.second)
        if IsDiscrete
            # first, discount proportional part to the floor-integer
            TmpVal = CF.second / ( 1.0 + IntRate * ( CF.first - floor(CF.first) ) )
            # then, discount it to floor(DiscountTo) moment
            TmpVal /= ( 1.0 + IntRate ) ^ ( floor(CF.first) - floor(DiscountTo) )
            # finally, discount to proportional DiscountTo moment
            Res = TmpVal * ( 1.0 + IntRate * ( DiscountTo - floor(DiscountTo) ) )
            return Res::Float64
        else
            Res = CF.second * exp( IntRate * ( DiscountTo - CF.first ) )
            return Res::Float64
        end
    elseif CF.first == Inf & isfinite(CF.second)
        return 0.0::Float64
    elseif ( CF.first == -Inf & CF.second == Inf ) | ( isfinite(CF.first) & CF.second == Inf )
        return Inf::Float64
    elseif ( isfinite(CF.first) & CF.second == -Inf ) | ( CF.first == -Inf & CF.second == -Inf )
        return -Inf::Float64
    else
        return NaN::Float64
    end
end # end functio PV
# ==================================
"""
    CashFlow

a synonym of Vector{Pair{Float64,Float64}},
whose element is (Time, Signed Payment);
denotes a series of money flows at different moments;
allows duplicated time/moments unless concise() it;

Mathematics:
1. The domain of Time is: {-Inf,R,Inf}, where R is finite real number
2. The domain of Money is: {-Inf,R,Inf}
3. then we define (mathematically) some operations on a single CashFlow
3. the principles of PV, please refer to the documentation of type Payment
"""
CashFlow = Vector{Payment}
# ----------------- Type Assertion
"""
    iscashflow( Obj::Any )

checks if an instance is one of CashFlow
"""
iscashflow( Obj::Any ) = isa( Obj, CashFlow )
# ----------------- Constructors
cashflow( Time::Real, Cash::Real ) = [Payment( Time, Cash )]
cashflow( Time::Real, Cash::Union{AbstractRange, Vector{T} where T <: Real} ) = begin
    local Res = CashFlow()  # initilize an instance to return
    for x in Cash; push!( Res, Payment( Time, x ) ); end;  # push new (time, payment) pairs to Res::CashFlow
    return Res::CashFlow
end
cashflow( Time::Union{AbstractRange, Vector{T} where T <: Real}, Cash::Real ) = begin
    local Res = CashFlow()
    for x in Time; push!( Res, Payment( x, Cash ) ); end;
    return Res::CashFlow
end
"""
    cashflow( Time::Vector{T} where T <: Real, Payment::Vector{T} where T <: Real )

creates a CashFlow from two vectors: one for time/moments, one for payments;
returns a CashFlow instance;
duplicated time/moments allowed;
"""
cashflow( Time::Union{AbstractRange, Vector{T} where T <: Real}, Cash::Union{AbstractRange, Vector{T} where T <: Real} ) = begin
    # valid & declare
    local T = length(Time)  # size of the cash flow
    T == length(Cash)  ||  error("expects two equal-length vectors but received: ( $(T), $(length(Cash)) )")
    local Res = CashFlow()  # initilize an instance to return
    for t in 1:T
        push!( Res, Payment( Time[t], Cash[t] ) )
    end
    return Res::CashFlow
end  # end function cashflow()
"""
    cashflow( TimeCash::Matrix{T} where T <: Real )

creates a CashFlow from a T*2 matrix (synonym of Array{T,2}),
where the first column is time/moments of each payment,
and the second column is payment amounts;
returns a CashFlow;
"""
cashflow( TimeCash::Matrix{T} where T <: Real ) = begin
    # valid & declare
    local T = size(TimeCash)  # size of the given matrix
    T[2] == 2  ||  error("expects a T*2 matrix but received: $(T)")
    local Res = CashFlow()
    for t in 1:T[1]; push!( Res, Payment( TimeCash[t,1], TimeCash[t,2] ) ); end;
    return Res::CashFlow
end  # end function cashflow()
"""
    cashflow( Time2Cash::Dict{T,T} where T <: Real )

creates a CashFlow from a number-to-number Dict (Hash table),
where the keys are time/moments of each payment,
and the values are payment amounts;
returns a CashFlow;
"""
cashflow( Time2Cash::Dict{T,T} where T <: Real ) = begin
    # valid & declare
    local Res = CashFlow()
    for x in Time2Cash; push!( Res, Payment( x.first, x.second ) ); end;
    return Res::CashFlow
end  # end function cashflow()
# ----------------- Operations
"""
    PV( CF::CashFlow, IntRate::Float64 ; DiscountTo::Union{Int,Float64} = 0.0, IsDiscrete::Bool = true )


computes the present value of a CashFlow at time CountTo::Union{Int,Float64} at a constant non-negative interest rate IntRate::Float64;

Receives:
1. CF::CashFlow, cash flow to discount
2. IntRate::Float64, a constant interest rate to discount, >=0, in digit rather than percentage
3. DiscountTo::Union{Int,Float64}, the moment to discount to
4. IsDiscrete::Bool, to discount as a discrete cash flow or a continuous cash flow, true by default
Returns:
1. Res::Float64, PV of the given cash flow
Depends on:
1. PV()

Note:
1. when IsDiscrete is false, money is discounted as `` PV = \\exp( - r (t-t_0) ) ``
1. when IsDiscrete is true, inner-year part is discounted at linear interpolation, e.g.
`` \\text{PV}_{t=1} = \\text{FV}_{t=1.23} \\cdot \\frac{1}{1+ 10% \\cdot 0.23 } ``
"""
function PV( CF::CashFlow, IntRate::Float64 ; DiscountTo::Real = 0.0, IsDiscrete::Bool = true )
    # NOTE: validations are did inside PV(::Payment,...)
    local Res::Float64 = 0.0  # result
    # compute
    for tmppayment in CF
        Res += PV( tmppayment, IntRate, DiscountTo = DiscountTo, IsDiscrete = IsDiscrete )
    end
    return Res::Float64
end # end function PV
# ----------------- Arithmetical Operations
# NOTE: + & - on CashFlow are set-level operations, not on element (Payment)
# NOTE: but - (get negative sign) works on single Payment
# NOTE: * & / are defined at element-level, not on set-level (CashFlow), please use .* and ./
"""
    Base.:+( cfa::CashFlow, cfb::CashFlow )

overloads addition (+) for CashFlow;
it is a kind of set operation (not defined on element!);
two CashFlow are unioned, where payments with different payment amounts but at the same moments are all kept;
and only one kept for those the same (time = time, amount = amount) Payments;
returns a new CashFlow instance;
"""
Base.:+( cfa::CashFlow, cfb::CashFlow ) = union( cfa, cfb )
"""
    Base.:-( cfa::CashFlow, cfb::CashFlow )

overloads minus (-) for CashFlows;
it is a kind of set operation (not defined on element!);
payments in cfb has opposite signs then added to cfa;
returns a new CashFlow instance;
"""
Base.:-( cfa::CashFlow, cfb::CashFlow ) = begin
    return ( cfa + ( .- cfb ) )::CashFlow
end  # end function -
# ----------------- Collection Operations
"""
    Base.cfcompress( cf::CashFlow )

compresses a cash flow, summing all duplicated records (same time but different payment amount),
and dropping zero payments (payment == 0);
returns a CashFlow with unique members;
"""
function cfcompress( cf::CashFlow )
    # declare
    local Res::CashFlow = CashFlow()
    #


end # end function unique



# gettimes, getpayments
# isfinite(::Payment)











# ==================================
"""
    AbstractAnnuity

a super abstract type for different annuities
"""
abstract type AbstractAnnuity <: Any end
# --------------------
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
struct StdAnnuity <: AbstractAnnuity
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
