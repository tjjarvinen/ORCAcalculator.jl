
# New AtomsCalculators style with state (execution) and parameters (medhod)
struct ORCAcalculatorbundle
    exection::OrcaExecutable
    method::OrcaMethod
end


AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(
    system,
    orca::ORCAcalculatorbundle;
    orca_stdout=stdout,
    ghosts=(),
    kwargs...
)
    res = calculate(system, orca.exection, orca.method; orca_stdout=orca_stdout, ghosts=ghosts)
    return res[:energy]
end


AtomsCalculators.@generate_interface function AtomsCalculators.forces(
    system, 
    orca::ORCAcalculatorbundle; 
    orca_stdout=stdout, 
    numgrad=false, 
    kwargs...
)
    res = calculate(system, orca.exection, orca.method; orca_stdout=orca_stdout, engrad=true, numgrad=numgrad)
    return res[:forces]
end


AtomsCalculators.promote_force_type(::Any, ::ORCAcalculatorbundle) = SVector(1., 1., 1.) * u"hartree/bohr" |> typeof