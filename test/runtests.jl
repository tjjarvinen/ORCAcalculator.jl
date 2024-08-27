using ORCAcalculator
using AtomsBase
using AtomsCalculators
using AtomsCalculators.Testing
using Unitful
using UnitfulAtomic
using Test

@testset "ORCAcalculator.jl" begin
    # Write your tests here.
    ox = ORCAexecutable()
    om = ORCAmethod("blyp def2-svp")
    orca = ORCAcalculatorbundle( ox, om )

    hydrogen = isolated_system([
        :H => [0, 0, 0.]u"Å",
        :H => [0, 0, 1.]u"Å"
    ])

    test_energy_forces(hydrogen, orca; orca_stdout=devnull)
    test_forces(hydrogen, orca; orca_stdout=devnull, numgrad=true)

    e_f = AtomsCalculators.energy_forces(hydrogen, orca)
    res = ORCAcalculator.get_results(orca)
    @test haskey(res, :dipolemoment)


    om = ORCAmethod("blyp def2-svvp") # Broken
    orca = ORCAcalculatorbundle( ox, om )
    @test_throws ProcessFailedException AtomsCalculators.potential_energy(hydrogen, orca; orca_stdout=devnull)
end
