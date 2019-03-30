"""
    EasyIO

user-friendly I/O methods. It covers:
1. timeTag(), generating time taggings
2. readcsv() & writecsv(), recoverage of two deprecated functions since Julia 0.6



"""
module EasyIO
   import Dates  # std lib, for time taggings
   import DelimitedFiles # std lib, for delimiter-separated file I/O

# -------------
"""
   timeTag()

generates a time tagging; the tag can be used to mark time or name files.
returns a String.
"""
function timeTag()
   local tagstr = replace(string(Dates.now()), "-" => "_" )
   tagstr = replace( tagstr, "T" => "_" )
   tagstr = replace( tagstr, ":" => "_" )
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




























end # end module EasyIO
