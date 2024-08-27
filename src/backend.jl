
"""
    ORCAexecutable(kwords...)

Holds information about ORCA executables and where calculation is performed.

Normally tries to find ORCA binary automatically. If this does not work, or you want to use
a special version. You can give path to the binary directly with `executable` keyword.

Different calculators can run on same directory when they have different `base_name`.

# Fields

- `executable::String`           :  path to orca executable binary
- `ncore::UInt=1`                :  number of cores used in calculations
- `memcore::UInt=1000`           :  maximum memory per cores in calcualations
- `tmp_dir::String=mktempdir()`  :  directory where calculations are performed
- `base_name::String`            :  base name for files that are created during calculation

# Example

Find orca executable and create temporary directory for calculations.
Calculations use two cores and max 2GB memory

```julia
ORCAexecutable(ncore=2, memcore=2000)
```

"""
struct ORCAexecutable{Int}
    "path for orca excecutable"
    executable::String
    "number of cores in calculation"
    ncore::UInt
    "maximum memory per core in megabytes"
    memcore::UInt
    "directory where calculations are performed"
    tmp_dir::String
    "prefix of files used in calculation"
    base_name::String
    function ORCAexecutable(;
        executable=nothing,
        ncore=1,
        maxmem=1000,
        tmp_dir=mktempdir(),
        base_name="orca_calculation",
        version=nothing
    )   
        
        if isnothing(version) && isnothing(executable)
            # look for orca_scf/orca_casscf instead of orca as there is an other Unix orca program
            # we can also determine version based on which binary is present
            try
                executable = readchomp(`which orca_scf`)[begin:end-4]
                version = 5
            catch _
                
            end
            if isnothing(version)
                try
                    executable = readchomp(`which orca_casscf`)[begin:end-7]
                    version = 6
                catch _
                    # we failed to find ORCA binary
                    error("Could not find ORCA binary")
                end
            end
            if ! isfile(executable)
                error("ORCA executable location does not have all ORCA binaries")
            end
        elseif isnothing(version) && isfile(executable) # Find ORCA version
            if isfile(executable*"_scf")
                version = 5
            elseif isfile(executable*"_casscf")
                version = 6
            else
                error("Could not determine ORCA version")
            end
        else
            error("Could not determine ORCA binary location or version")
        end
        if ! isdir(tmp_dir)
            error("tmp dir does not exist")
        end
        if isnothing(version)
            error("ORCA version is undefined")
        end
        new{version}(executable, ncore, maxmem, tmp_dir, base_name)
    end
end


"""
    ORCAmethod(method::AbstractString, control::AbstractString="")

Used to write method definition to ORCA input file.

`method` will be prefixed with "!" in input file.
This is mainly ment to define calculation method and basis.

`control` is added as is. This allows you to put in general input
like control blocks starting with "%"

# Example

A simple DFT calculation

```julia
ORCAmethod("B3LYP D4 def2-TZVPP")
```

Use control block to control calculation

```julia
ORCAmethod(
    "RI-MP2 SVP def2-SVP/C",
    "%mp2 natorbs true
        density unrelaxed
    end"
)
```
"""
struct ORCAmethod
    "will be prefixed with ! in ORCA input"
    method::String
    "will be added as is to ORCA input"
    control::String
    function ORCAmethod(method::AbstractString, control::AbstractString="")
        new(method, control)
    end
end


function write_input(
    io, 
    system, 
    oex::ORCAexecutable, 
    oct::ORCAmethod; 
    add_engrad=false, 
    add_numgrad=false,
    ghosts=()
)
    if add_engrad
        println(io, "! ENGRAD")
    end
    if add_numgrad
        println(io, "! NUMGRAD")
    end
    println(io, "! ", oct.method)
    println(io, oct.control)
    println(io, "! MINIPRINT\n")
    if oex.ncore > 1
        println(io, "%pal nprocs $(oex.ncore) end")
    end
    println(io, "%maxcore $(oex.memcore)\n")
    println(io, "*xyz 0 1")
    foreach( enumerate(system) ) do (i,atom)
        s = atomic_symbol(atom)
        r = ustrip.( u"Å", position(atom) )
        if i in ghosts
            println(io, s, " : ", r[1], " ", r[2], " ", r[3])
        else
            println(io, s, " ", r[1], " ", r[2], " ", r[3])
        end
    end
    println(io, "*")
end

"""
    calculate(...)

Does the actual calculations

"""
function calculate(
    system, 
    oex::ORCAexecutable, 
    oct::ORCAmethod; 
    orca_stdout=devnull, 
    engrad=false, 
    numgrad=false, 
    ghosts=(),
    clean_files=true
)
    # Clean old files
    if clean_files
        files = readdir(oex.tmp_dir)
        foreach(files) do file
            if occursin( Regex(oex.base_name * "*"), file )
                rm( joinpath(oex.tmp_dir, file) )
            end
        end
    end

    # Write input
    f_input = joinpath(oex.tmp_dir, oex.base_name * ".inp")
    open(f_input, "w") do io
        write_input(io, system, oex, oct; add_engrad=engrad, add_numgrad=numgrad, ghosts=ghosts)
    end

    # Perform calculation
    try
        (run ∘ pipeline)(`$(oex.executable) $f_input`; stdout=orca_stdout)
    catch err
        # Remove temporary files that can be very large
        foreach(readdir(oex.tmp_dir)) do file
            if occursin( Regex(oex.base_name * "*.tmp"), file )
                rm( joinpath(oex.tmp_dir, file) )
            end
        end
        throw(err)
    end

    return nothing
end


function get_results(oex::ORCAexecutable)
    tmp1 = Threads.@spawn parse_engrad_file(oex)
    tmp2 = Threads.@spawn parse_property_file(oex)
    d1 = fetch(tmp1)
    d2 = fetch(tmp2)
    foreach(pairs(d2)) do (k,v)
        d1[k] = v
    end
    return d1
end



function parse_engrad_file(oex::ORCAexecutable)
    orca_engrad_file = joinpath(oex.tmp_dir, oex.base_name * ".engrad")
    out = Dict{Symbol, Any}()
    if isfile(orca_engrad_file)
        res = read(orca_engrad_file, String)
        lines = split(res, '\n')

        natoms = parse(Int, lines[4])
        
        e = parse(Float64, lines[8])
        out[:energy] = e * hartree

        f_vector = map( 12:12+3natoms-1 ) do i
            parse(Float64, lines[i])
        end
        f = reinterpret(SVector{3,Float64}, f_vector) .* (hartree/bohr)
        out[:forces] = f
    end
    return out
end


function parse_property_file(oex::ORCAexecutable{6})
    out = Dict{Symbol, Any}()
    orca_results_file = joinpath(oex.tmp_dir, oex.base_name * ".property.txt")
    if isfile(orca_results_file)
        res = read(orca_results_file, String)
        lines = split(res, "\n")
        
        # Energy
        rgx_energy = r"TOTALENERGY\s+.*\s+(?<energy>[+\-]?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+\-]?\d+)?)\s+\"Hartrees\""
        for line in lines
            if occursin(rgx_energy, line)
                m = match(rgx_energy, line)
                e = parse(Float64, m[:energy])
                out[:energy] = e * hartree
                break
            end
        end

        # Dipolement
        rgx_dipole =  r"&DIPOLETOTAL \[&Type \"ArrayOfDoubles\", &Dim \(3,1\)\] \"Total\""
        dipole_lines = nothing
        for (i,line) in enumerate(lines)
            if occursin(rgx_dipole, line)
                dipole_lines = i+3:i+5
            end
        end
        if !isnothing(dipole_lines)
            dipoles = [ parse(Float64, split(lines[i])[2]) for i in dipole_lines  ]
            out[:dipolemoment] = SVector(dipoles...) * debye
        end
    end
    return out
end


function parse_property_file(oex::ORCAexecutable{5})
    out = Dict{Symbol, Any}()
    orca_results_file = joinpath(oex.tmp_dir, oex.base_name * "_property.txt")
    if isfile(orca_results_file)
        res = read(orca_results_file, String)
        lines = split(res, "\n")

        # Energy
        for line in lines
            if occursin(r"Total Energy*\d*", line)
                e = parse(Float64, split(line)[3])
                out[:energy] = e * hartree
                break
            end
        end

        # Dipolement
        dipole_lines = nothing
        for (i,line) in enumerate(lines)
            if occursin(r"Total Dipole moment:", line)
                dipole_lines = i+2:i+4
            end
        end
        dipoles = [ parse(Float64, split(lines[i])[2]) for i in dipole_lines  ]
        out[:dipolemoment] = SVector(dipoles...) .* debye
    end
    return out
end