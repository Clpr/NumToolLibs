"""
    EasyIO

user-friendly I/O methods. It covers:
1. timeTag(), generating time taggings
2. readcsv() & writecsv(), recoverage of two deprecated functions since Julia 0.6



"""
module EasyIO
   import LinearAlgebra: diagind, diag  # std lib, diagind() gets index range for an assigned offline parameter, used in diagwrite()
   import Dates  # std lib, for time taggings
   import DelimitedFiles # std lib, for delimiter-separated file I/O
   export timeTag, readcsv, writecsv, diagw, diagw!, expand, expand!

# -------------
"""
   timeTag()

generates a time tagging; the tag can be used to mark time or name files.
returns a String.
"""
function timeTag()
   local tagstr = replace(string(Dates.now()), ":" => "_" )
   # tagstr = replace( tagstr, "T" => "_" )
   tagstr = replace( tagstr, "." => "_" )
   return tagstr::String
end
# -------------
# Two IO functions from Julia 0.6
# NOTE: well ... the readcsv() & writecsv() before Julai 0.6 are SO convenient in some ways ...
#       I re-define the two functions here; (directly copied from base Julia 0.6)
readcsv(io; opts...)          = DelimitedFiles.readdlm(io, ','; opts...)
readcsv(io, T::Type; opts...) = DelimitedFiles.readdlm(io, ',', T; opts...)
"""
    writecsv(filename, A; opts)
Equivalent to [`writedlm`](@ref) with `delim` set to comma.
"""
writecsv(io, a; opts...) = DelimitedFiles.writedlm(io, a, ','; opts...)
# -------------
"""
    diagw( MAT::Matrix, VEC::Union{Vector,Tuple,NamedTuple,AbstractRange} ; offset::Int = 0 )

Writes a vector "VEC" into a matrix "MAT" diagonally, according to offset= (value starts from 0).
Returns a copy of MAT;

Depends:
1. LinearAlgebra.diagind()

Note:
1. diagind(MAT) (the index of positions to fill) can have different length from vector VEC
   1. if length(diagind(MAT)) < length(VEC): extra elements of VEC will be neglected
   2. if length(diagind(MAT)) > length(VEC): extra elements of MAT[diagind(MAT)] will be kept
"""
function diagw( MAT::Matrix, VEC::Union{Vector,Tuple,NamedTuple,AbstractRange} ; offset::Int = 0 )
   # declare
   local Res::Matrix = copy(MAT)
   local tmpIdxRange = diagind(Res, offset) # get indexes range
   # fill elements
   for s in 1:min( length(VEC), length(tmpIdxRange) ); Res[ tmpIdxRange[s] ] = VEC[s]; end
   return Res::Matrix
end # end function diagw
# ----------------
"""
    diagw!( MAT::Matrix, VEC::Union{Vector,Tuple,NamedTuple,AbstractRange} ; offset::Int = 0 )

Writes a vector "VEC" into a matrix "MAT" diagonally (in-place), according to offset= (value starts from 0).
Returns nothing;

Depends:
1. LinearAlgebra.diagind()

Note:
1. diagind(MAT) (the index of positions to fill) can have different length from vector VEC
   1. if length(diagind(MAT)) < length(VEC): extra elements of VEC will be neglected
   2. if length(diagind(MAT)) > length(VEC): extra elements of MAT[diagind(MAT)] will be kept
"""
function diagw!( MAT::Matrix, VEC::Union{Vector,Tuple,NamedTuple,AbstractRange} ; offset::Int = 0 )
   # declare
   local tmpIdxRange = diagind(MAT, offset) # get indexes range
   # fill elements
   for s in 1:min( length(VEC), length(tmpIdxRange) ); MAT[ tmpIdxRange[s] ] = VEC[s]; end
   return nothing::Nothing
end # end function diagw
# ----------------
"""
   expand( VEC::Vector, LEN::Int )

expands VEC to specific length LEN;
if LEN > length(VEC), fill extra positions with the last element of VEC;
if LEN < length(VEC), drop the ending elements;
returns a Vector/1-D Array of the same type
"""
function expand( VEC::Vector, LEN::Int )
   LEN >= 0  ||  "expects a non-negative length to expand but received: $(LEN)"
   tmpLen = length(VEC); tmpVec = copy(VEC)
   if LEN > tmpLen
      tmpVal = tmpVec[end]
      for s in range( tmpLen+1, length= LEN-tmpLen )
         push!( tmpVec, tmpVal )
      end
      return tmpVec::Vector
   else
      return tmpVec[1:LEN]::Vector
   end
end # end function expand
function expand!( VEC::Vector, LEN::Int )
   LEN >= 0  ||  "expects a non-negative length to expand but received: $(LEN)"
   tmpLen = length(VEC)
   if LEN > tmpLen
      tmpVal = VEC[end]
      for s in range( tmpLen + 1, length= LEN - tmpLen )
         push!( VEC, tmpVal )
      end
      return nothing::Nothing
   else
      splice!( VEC, (LEN + 1):tmpLen ) # similar to pop!() but on a range
      return nothing::Nothing
   end
end # end function expand!
# ----------------

























end # end module EasyIO
