# This contains a few scripts to run Perple_X examples


"""
    success = run_test_pipe(; dir="bl691_MEEMUM", inputfile="bl691_MEEMUM_input.txt", dir_exec="../sources/", meemum=false, vertex=false)  

Runs a test case in directory `dir`, using input file `inputfile` that is piped into either the `vertex` or `meemum` executable. 
The executable is located in directory `dir_exec` relative to the current directory.

If something goes wrong, it will return false.
"""
function run_test_pipe(; dir="bl691_MEEMUM", inputfile="bl691_MEEMUM_input.txt", dir_exec="../src/", meemum=false, vertex=false)  
    cur_dir = pwd();

    # Indicate the executable you will run
    exec=nothing
    if meemum==true
        exec = "meemum"
    elseif vertex==true
        exec = "vertex"
    end
    if isnothing(exec)
        error("You need to indicate which executable you want to run")
    end

    # Try running the executable while piping the input file; in case of error it'll return false
    success=false
    try 
        dir_pipe = joinpath(cur_dir,dir)
        cd(dir) # change to the directory
        p = pipeline(`../$dir_exec/$exec`, stdin=inputfile)       # pipe input file to executable
        @show p
        run(p)  # run
        success=true        
    catch
        success=false   # something errored
    end
    cd(cur_dir)
    return success
end




"""

Reads a meemum output file 
"""
function read_meemum(;filename="bl691_MEEMUM", dir="bl691_MEEMUM")

    cur_dir = pwd();
    file    = joinpath(cur_dir,dir,filename*".prn")

    # Read floats from file
    T,P, 
    Cp, spec_Cp,
    S,
    H = extract_value(file, ("T(K)", "P(bar)",
                                            "Heat Capacity (J/K/kg)","Specific Heat Capacity (J/K/m3)",
                                            "Entropy (J/K/kg)",
                                            "Enthalpy (J/kg)"
                                            )) 

    variance = extract_value(file, "Variance (c-p+2)")

    seismic_props  = extract_data_block(file,"Seismic Properties:")
    mol_props      = extract_data_block(file,"Molar Properties and Density:")
    phase_comp     = extract_data_block(file,"Phase Compositions (molar  proportions):")
    chem_pot       = read_chemical_potentials(file);
    
    cd(cur_dir)

    # put all data in a Named Tuple
    out = (T=T, P=P, Cp=Cp, spec_Cp=spec_Cp, S=S, H=H, variance=variance, 
                chem_pot=chem_pot, seismic_props=seismic_props, phase_comp=phase_comp)

    return out
end


"""
    extract_value(file::String, keywords::NTuple{N,String}) 

This scans the file `file` for `keywords` and returns the numerical extract_value 
"""
function  extract_value(file::String, keyword::NTuple{N,String}) where N
    
    d =  Vector{Union{Nothing, Float64}}(nothing, N)
    open(file) do f
        while ! eof(f) 
            # read a new / next line for every iteration          
            line = readline(f)
            for i=1:N        
                if contains(line, keyword[i])
                    d[i] = parse(Float64,split(line,'=')[2]); 
                end
            end
        end
    end

    return Tuple(d)
end


"""
    extract_value(file::String, keyword::String)

Extract a single value of a line
"""
function  extract_value(file::String, keyword::String)
    return extract_value(file, (keyword,))[1]
end


"""

    This extracts a data block from a figure
"""
function  extract_data_block(file::String, keyword::String, additional_empty_lines=0)
    data = [];
    open(file) do f
        while ! eof(f) 
            line = readline(f)
            if contains(line, keyword)

                line_names = readline(f)
                line_names = replace(line_names, "Poisson ratio"=>"Poisson_ratio")

                names = split(remove_space_beforeafter_symbol(line_names,"%"));
                sym_names = Tuple([:Name; Symbol.(names)]);
      
                # read the block itself
                while !isempty(line)
                    line = readline(f)
                    for i=1:additional_empty_lines
                        line = readline(f)
                    end
                    line = replace(line, " - "=>"_")
                    
                    if !isempty(line)

                        # convert values to float
                        values = [split(line)[1]; parse.(Float64,split(line)[2:end])];
                        

                        # Store as NamedTuple
                        data_line = NamedTuple{sym_names}(Tuple(values))

                        push!(data, data_line)
                    end
                end
          
            end
        end
    end

    return data
end
    
function remove_space_beforeafter_symbol(str::String, sym::String)
    
    ind_all = findall(sym, str)
    id  = [ind_all[i][1]-1 for i=1:length(ind_all)];
    id1 = [ind_all[i][1]+1 for i=1:length(ind_all)];
    
    id_vec = setdiff(Vector(1:length(str)), id)
    id_vec = setdiff(id_vec, id1)

    str = str[id_vec]
    return str
end


"""
    out = read_chemical_potentials(file::String)
reads chemical potentials from a meemum file
"""
function  read_chemical_potentials(file::String)
    
    out = nothing
    open(file) do f
        while ! eof(f) 
            # read a new / next line for every iteration          
            line = readline(f)    
            if contains(line, "Chemical Potentials (J/mol):")
                line = readline(f)  
                names = Symbol.(split(readline(f)))
                values = parse.(Float64,split(readline(f)))

                out = NamedTuple{Tuple(names)}(Tuple(values))
            end

        end
    end

    return out
end
