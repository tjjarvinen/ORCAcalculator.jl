module ORCAcalculator

using AtomsBase
using AtomsCalculators
using Compat
using StaticArrays
using Unitful
using UnitfulAtomic

export ORCAexecutable
export ORCAmethod
export ORCAcalculatorbundle

@compat public get_results
@compat public calculate


# Define dipolomement unit to ease use
@unit debye "D" debye (1e-21/299792458)*u"C*m" false

# Define context untits for better conversions
const bohr = Unitful.ContextUnits(u"bohr", u"Å")
const hartree = Unitful.ContextUnits(u"hartree", u"eV")

function __init__()
    Unitful.register(@__MODULE__)
end

include("backend.jl")
include("atoms_calculators.jl")



end # module
