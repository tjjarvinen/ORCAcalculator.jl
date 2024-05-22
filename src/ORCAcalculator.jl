module ORCAcalculator

using AtomsBase
using AtomsCalculators
using StaticArrays
using Unitful
using UnitfulAtomic

export ORCAexecutable
export ORCAmethod
export ORCAcalculatorbundle

@unit debye "D" debye (1e-21/299792458)*u"C*m" false

function __init__()
    Unitful.register(@__MODULE__)
end

include("backend.jl")
include("atoms_calculators.jl")



end # module
