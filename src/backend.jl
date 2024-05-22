

struct ORCAexecutable
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
    function ORCAexecutable(;
        executable=nothing,
        ncore=1,
        maxmem=1000,
        tmp_dir=mktempdir(),
        base_name="orca_calculation"
    )   
        orca_location = something(executable, readchomp(`which orca`))
        if ! isfile(orca_location*"_scf")
            error("ORCA executable location does not have all ORCA binaries")
        end
        if ! isdir(tmp_dir)
            error("tmp dir does not exist")
        end
        new(orca_location, ncore, maxmem, tmp_dir, base_name)
    end
end

struct ORCAmethod
    method::String
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


function calculate(system, oex::ORCAexecutable, oct::ORCAmethod; orca_stdout=stdout, engrad=false, numgrad=false, ghosts=())
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
        write_input(io, system, oex, oct; add_engrad=engrad, add_numgrad=numgrad, ghosts=ghosts)
    end

    # Perform calculation
    try
        (run ∘ pipeline)(`$(oex.executable) $f_input`; stdout=orca_stdout)
    catch err
        foreach(files) do file
            if occursin( Regex(oex.base_name * "*.tmp"), file )
                rm( joinpath(oex.tmp_dir, file) )
            end
        end
        throw(err)
    end

    # Read results
    out = Dict{Symbol, Any}()
    orca_results_file = joinpath(oex.tmp_dir, oex.base_name * "_property.txt")
    if isfile(orca_results_file)
        res = read(orca_results_file, String)

        # Energy
        for line in eachline(IOBuffer(res))
            if occursin(r"Total Energy*\d*", line)
                e = parse(Float64, split(line)[3])
                out[:energy] = e * u"hartree"
                break
            end
        end

        # Dipolement
        lines = split(res, "\n")
        dipole_lines = nothing
        for (i,line) in enumerate(lines)
            if occursin(r"Total Dipole moment:", line)
                dipole_lines = i+2:i+4
            end
        end
        dipoles = [ parse(Float64, split(lines[i])[2]) for i in dipole_lines  ]
        out[:dipolemoment] = SVector(dipoles...) .* debye
    end

    ## Forces output ##
    ###################
    orca_engrad_file = joinpath(oex.tmp_dir, oex.base_name * ".engrad")
    if isfile(orca_engrad_file)  # we have .engrad file

        res = read(orca_engrad_file, String)

        # We need to look for the block type
        # 0 is unknown block
        # 1 is energy block
        # 2 is gradient block coming next, but numbers not started yet
        # 3 is gradient block numbers started
        t = 0
        grad = Float64[]
        sizehint!(grad, 3*length(system))
        for line in eachline(IOBuffer(res))
            #@info "t = $t"
            if t == 0  # look for block type
                if occursin(r"The current gradient in Eh/bohr", line)
                    t = 2
                elseif occursin(r"The current total energy in Eh", line)
                    t = 1
                end
            else
                if t == 1 && line[begin] != '#' # energy line number
                    out[:energy] = parse(Float64, line) * u"hartree"
                    t = 0
                elseif t == 2 && line[begin] != '#' # gradient line first number
                    f = parse(Float64, line)
                    append!(grad, f)
                    t = 3
                elseif t == 3 && line[begin] != '#' # gradient line numbers from 2->
                    f = parse(Float64, line)
                    append!(grad, f)
                elseif t == 3 && line[begin] == '#'
                    #all done
                    break
                end
            end
        end
        tmp = reshape(grad, 3, :)
        out[:forces] = - reinterpret(reshape, SVector{3, Float64}, tmp) * u"hartree/bohr"
    end
    return out
end