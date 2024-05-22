
# New AtomsCalculators style with state (execution) and parameters (medhod)

"""
    ORCAcalculatorbundle

Used to wrap `ORCAexecutable` and `ORCAmethod` together to form a AtomsCalculators
compatable calculator.

# Fields
- exection::ORCAexecutable
- method::ORCAmethod

# Example

Simple calculation using AtomsCalculators

```julia
using ORCAcalculator
using AtomsBase
using AtomsCalculators
using Unitful

orca = ORCAcalculatorbundle(
    ORCAexecutable()
    ORCAmethod("blyp def2-svp")
)
```
hydrogen = isolated_system([
    :H => [0, 0, 0.]u"Å",
    :H => [0, 0, 1.]u"Å"
])

AtomsCalculators.potential_energy(hydrogen, orca)
"""
struct ORCAcalculatorbundle
    exection::ORCAexecutable
    method::ORCAmethod
end


AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
    ::AtomsCalculators.Energy,
    system,
    orca::ORCAcalculatorbundle;
    orca_stdout=devnull,
    ghosts=(),
    kwargs...
)
    res = calculate(system, orca.exection, orca.method; orca_stdout=orca_stdout, ghosts=ghosts)
    return NamedTuple( pairs(res) )
end


AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
    ::AtomsCalculators.Forces,
    system, 
    orca::ORCAcalculatorbundle; 
    orca_stdout=devnull, 
    numgrad=false, 
    kwargs...
)
    res = calculate(system, orca.exection, orca.method; orca_stdout=orca_stdout, engrad=true, numgrad=numgrad)
    return NamedTuple( pairs(res) )
end


AtomsCalculators.promote_force_type(::Any, ::ORCAcalculatorbundle) = SVector(1., 1., 1.) * u"hartree/bohr" |> typeof


function AtomsCalculators.energy_forces(
    system, 
    orca::ORCAcalculatorbundle; 
    orca_stdout=devnull, 
    numgrad=false, 
    kwargs...
)
    res = calculate(system, orca.exection, orca.method; orca_stdout=orca_stdout, engrad=true, numgrad=numgrad)
    return NamedTuple( pairs(res) )
end