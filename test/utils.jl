# This contains a few scripts to run Perple_X examples

"""
success = run_meemum(; dir="bl691_MEEMUM", inputfile="bl691_MEEMUM_input.txt", dir_meemum="../sources/meemum")  

Runs a meemum test case in directory `dir`, using input file `inputfile`. The `meemum` executable is located in directory `dir_meemum` relative to the current directory.

If something goes wrong, it will return false
"""
function run_meemum(; dir="bl691_MEEMUM", inputfile="bl691_MEEMUM_input.txt", dir_meemum="../sources/meemum")  
cur_dir = pwd();

# Try running; in case of error it'll return false
success=false
try 
    cd(dir)
    p = pipeline(`../$dir_meemum`, stdin=inputfile)
    run(p)
    success=true
catch
    success=false
end
cd(cur_dir)
return success
end


