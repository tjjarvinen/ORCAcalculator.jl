using ORCAcalculator
using AtomsBase
using AtomsCalculators
using AtomsCalculators.AtomsCalculatorsTesting
using Unitful
using UnitfulAtomic
using Test

@testset "ORCAcalculator.jl" begin
    # Write your tests here.
    ox = OrcaExecutable()
    om = OrcaMethod("blyp", "def2-svp")
    orca = ORCAcalculatorbundle( ox, om )

    hydrogen = isolated_system([
        :H => [0, 0, 0.]u"Å",
        :H => [0, 0, 1.]u"Å"
    ])

    test_energy_forces(hydrogen, orca; orca_stdout=devnull)
end
