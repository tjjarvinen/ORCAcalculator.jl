module ORCAcalculator

using AtomsBase
using AtomsCalculators
using Unitful
using UnitfulAtomic

# Write your package code here.

const supported_output = [
    :energy,
    :forces,
    :num_forces,
    :dipolemoment
]

struct OrcaExecutable
    "path for orca excecutable"
    executable
    "number of cores in calculation"
    ncore::UInt
    "maximum memory per core in megabytes"
    memcore::UInt
    "directory where calculations are performed"
    tmp_dir
    "prefix of files used in calculation"
    base_name
    function OrcaExecutable(;
        executable=nothing,
        ncore=1,
        maxmem=1000,
        tmp_dir=mktempdir(),
        base_name="orca_calculation"
    )   
        orca_location = something(executable, readchomp(`which orca`))
        if ! isdir(tmp_dir)
            error("tmp dir does not exist")
        end
        new(orca_location, ncore, maxmem, tmp_dir, base_name)
    end
end

struct OrcaCalculationType
    basis::String
    method::String
    orca_results::Vector{Symbol}
    function OrcaCalculationType(basis::AbstractString, method::AbstractString, orca_results::Symbol...)
        foreach( orca_results ) do x
            if ! (x in supported_output)
                error("Not supported ORCA output: $x")
            end
        end
        new(basis, method, [x for x in orca_results])
    end
end


function write_input(io, system, oex::OrcaExecutable, oct::OrcaCalculationType)
    println(io, "! ", oct.method)
    println(io, "! ", oct.basis)
    println(io, "! MINIPRINT\n")
    if oex.ncore > 1
        println(io, "%pal nprocs $(oex.ncore) end")
    end
    println(io, "%maxcore $(oex.memcore)\n")
    println(io, "*xyz 0 1")
    foreach(system) do atom
        s = atomic_symbol(atom)
        r = ustrip.( u"Å", position(atom) )
        println(io, s, " ", r[1], " ", r[2], " ", r[3])
    end
    println(io, "*")
end


function calculate(system, oex::OrcaExecutable, oct::OrcaCalculationType; orca_stdout=stdout)
    # Clean old files
    files = readdir(oex.tmp_dir)
    foreach(files) do file
        if occursin( Regex(oex.base_name * "*"), file )
            rm( joinpath(oex.tmp_dir, file) )
        end
    end

    # Write input
    f_input = joinpath(oex.tmp_dir, oex.base_name * ".inp")
    open(f_input, "w") do io
        write_input(io, system, oex, oct)
    end

    # Perform calculation
    orca_status = (run ∘ pipeline)(`$(oex.executable) $f_input`; stdout=orca_stdout)
    if orca_status.termsignal != 0
        # delete temporary files and raise error
        foreach(files) do file
            if occursin( Regex(oex.base_name * "*.tmp"), file )
                rm( joinpath(oex.tmp_dir, file) )
            end
        end
        error("there was an problem with orca call")
    end

    # Read results
    orca_results_file = joinpath(oex.tmp_dir, oex.base_name * "_property.txt")
    res = read(orca_results_file, String)
    out = Dict{Symbol, Any}()
    for line in eachline(IOBuffer(res))
        if occursin(r"Total Energy*\d*", line)
            e = parse(Float64, split(line)[3])
            out[:energy] = e * u"hartree"
            break
        end
    end
    return out
end


end
