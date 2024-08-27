
# New AtomsCalculators style with state (execution) and parameters (method)

"""
    ORCAcalculatorbundle

Used to wrap `ORCAexecutable` and `ORCAmethod` together to form a AtomsCalculators
compatable calculator.

# Fields
- execution::ORCAexecutable
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

hydrogen = isolated_system([
    :H => [0, 0, 0.]u"Å",
    :H => [0, 0, 1.]u"Å"
])

AtomsCalculators.potential_energy(hydrogen, orca)
```

You can also use low level interface

```julia
AtomsCalculators.calculate(
    AtomsCalculators.Energy(),
    hydrogen,
    orca, 
    ORCAmethod(HF def2-svp), 
    ORCAexecutable()
)
```
"""
struct ORCAcalculatorbundle
    execution::ORCAexecutable
    method::ORCAmethod
end

get_results(orca::ORCAcalculatorbundle) = get_results(orca.execution)
parse_engrad_file(orca::ORCAcalculatorbundle) = parse_engrad_file(orca.execution)
parse_property_file(orca::ORCAcalculatorbundle) = parse_property_file(orca.execution)


AtomsCalculators.energy_unit(::ORCAcalculatorbundle) = hartree
AtomsCalculators.length_unit(::ORCAcalculatorbundle) = bohr

AtomsCalculators.get_state(oeb::ORCAcalculatorbundle) = oeb.execution
AtomsCalculators.get_parameters(oeb::ORCAcalculatorbundle) = oeb.method

AtomsCalculators.set_state!(oeb::ORCAcalculatorbundle, ox::ORCAexecutable) = ORCAcalculatorbundle(ox, oeb.method) 
AtomsCalculators.set_parameters!(oeb::ORCAcalculatorbundle, om::ORCAmethod) = ORCAcalculatorbundle(oeb.execution, om)


AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
    ::AtomsCalculators.Energy,
    system,
    orca::ORCAcalculatorbundle,
    pr::Union{Nothing,ORCAmethod}=orca.method,
    st::Union{Nothing,ORCAexecutable}=orca.execution;
    orca_stdout=devnull,
    ghosts=(),
    kwargs...
)
    if isnothing(pr)
        pr = orca.method
    end
    if isnothing(st)
        st = orca.execution
    end
    calculate(system, st, pr; orca_stdout=orca_stdout, ghosts=ghosts)
    res = parse_property_file(orca)
    res[:state] = AtomsCalculators.get_state(orca)
    return NamedTuple( pairs(res) )
end


AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
    ::AtomsCalculators.Forces,
    system, 
    orca::ORCAcalculatorbundle,
    pr::Union{Nothing,ORCAmethod}=orca.method,
    st::Union{Nothing,ORCAexecutable}=orca.execution; 
    orca_stdout=devnull, 
    numgrad=false, 
    kwargs...
)
    if isnothing(pr)
        pr = orca.method
    end
    if isnothing(st)
        st = orca.execution
    end
    calculate(system, st, pr; orca_stdout=orca_stdout, engrad=true, numgrad=numgrad)
    res = parse_engrad_file(orca)
    res[:state] = AtomsCalculators.get_state(orca)
    return NamedTuple( pairs(res) )
end

# We only have one output type for forces
AtomsCalculators.promote_force_type(::AtomsBase.AbstractSystem, ::ORCAcalculatorbundle) = SVector(1., 1., 1.) * (hartree/bohr) |> typeof


function AtomsCalculators.energy_forces(
    system, 
    orca::ORCAcalculatorbundle; 
    orca_stdout=devnull, 
    numgrad=false, 
    kwargs...
)
    calculate(system, orca.execution, orca.method; orca_stdout=orca_stdout, engrad=true, numgrad=numgrad)
    res = parse_engrad_file(orca)
    return NamedTuple( pairs(res) )
end