# ORCAcalculator

Julia inteface to [ORCA](https://www.faccts.de/orca/).

Currently supported calculation results are:

- Energy
- Forces with analytical and numerical gradients
- Dipolemoment

## Installation

You need to download ORCA from the official [forum](https://orcaforum.kofo.mpg.de) and install it.
If you add ORCA directory to `PATH` then everything should work.
The package will call `which orca` at creation of `ORCAexecutable` and, if the given path then includes
other ORCA executables, you should be good to go. If not, you can give path to ORCA executable
on creation of `ORCAexecutable`

```julia
ORCAexecutable( executable="path to orca" )
```

## Running calculations

To perform calculation you need to create two structures:
1. `ORCAexecutable` that has information of the binaries,
where calculations is performed and how many cores and how much memory is used.
2. `ORCAmethod` that defines what calculation method and basis set is used.

To create `ORCAexecutable` you can call

```julia
oex = ORCAexecutable(
    executable= "path to orca binary", # will be seached if not given
    ncore=1,          # number of cores used in calculations
    maxmem=1000,      # how much memory in MB is used per core
    tmp_dir=mktempdir(),  # directory where calculations are performed
    base_name="orca_calculation" # this controls file naming
)
```

To create `ORCAmethod` you can call

```julia
ORCAmethod(
    method::String,     # method like "blyp def2-svp"
    control::String=""  # additional control block - meant for block input
)
```

After this you can call `ORCAcalculator.calculate` to perform calculations on [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) structures

```julia
ORCAcalculator.calculate(system, orca)
```

## AtomsCalculators Support

[AtomsCalculators interface](https://github.com/JuliaMolSim/AtomsCalculators.jl) is sopported with `ORCAcalculatorbundle` structure that combines `ORCAexecutable` and `ORCAmethod`

```julia
ORCAcalculatorbundle(
    ORCAexecutable(),
    ORCAmethod("blyp def2-svp")
)
```

that can then be used an calculator in programs that suport AtomsCalculators interface.

### Example

```julia
using AtomsBase
using AtomsCalculators
using ORCAcalculator
using Unitful

hydrogen = isolated_system([
    :H => [0, 0, 0.]u"Å",
    :H => [0, 0, 1.]u"Å"
])

orca = ORCAcalculatorbundle(
    ORCAexecutable(),
    ORCAmethod("blyp def2-svp")
)

AtomsCalculators.potential_energy(hydrogen, orca)

# Print ORCA output
AtomsCalculators.potential_energy(hydrogen, orca; orca_stdout=stdout)


# Forces
AtomsCalculators.forces(hydrogen, orca)

# Forces using numerical gradient
AtomsCalculators.forces(hydrogen, orca; numgrad=true)
```

You can get access to dipolemoment by using calculator interface

```julia
# Works for Energy() and Forces()
res = AtomsCalculators.calculate( AtomsCalculators.Energy(), hydrogen, orca )

@show res[:dipolemoment]
```

Also combination calls return dipole moment.

### Low level interface

Low level (calculator) interface uses `ORCAmethod` as `parameter` and `ORCAexecutable` as `state`

```julia
res = AtomsCalculators.calculate(
    AtomsCalculators.Energy(),
    hydrogen,
    orca, 
    ORCAmethod(HF def2-svp), 
    ORCAexecutable()
)

res[:state] # this is ORCAexecutable
```

Currently `ORCAexecutable` does not hold any information about the system.
In the future new functionality to access e.g. orbitals could be added.

## Example calculate non-covalent interraction between two nitrogen molecules

This example uses counterpoise correction to
correct basis set superposition error between two molecules.

```julia
using AtomsBase
using AtomsCalculators
using ORCAcalculator
using Unitful


orca = ORCAcalculatorbundle(
    ORCAexecutable(ncore=4, maxmem=2000),
    ORCAmethod("CCSD(T)-F12/RI TIGHTSCF cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C")
)

n2_complex = isolated_system([
    :N => [0.0, 0.0, 0.0]u"Å",
    :N => [0.0, 0.0, 1.1]u"Å",

    :N => [3.5, 0.0, 0.0]u"Å",
    :N => [4.6, 0.0, 0.0]u"Å",
])


E12 = AtomsCalculators.potential_energy(n2_complex, orca)

# ghosts add basis functions to atoms 1 and 2, while only 3 and 4 have nucleus and electrons
E2 = AtomsCalculators.potential_energy(n2_complex, orca; ghosts=[1,2])

E1 = AtomsCalculators.potential_energy(n2_complex, orca; ghosts=[3,4])

E = E12 - E1 - E2
```
