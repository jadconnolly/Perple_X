using Test

include("utils.jl")

@testset "MEEMUM tests" begin

    dir_test = "bl691_MEEMUM";
    dir_expected = "bl691_MEEMUM/output";
    file_test = "bl691_MEEMUM"
    file_input = "bl691_MEEMUM_input.txt"

    out = run_test_pipe(dir=dir_test, inputfile=file_input, meemum=true)  
    @test out == true

    # Read expected & new files
    out_expected = read_meemum(dir=dir_expected, filename=file_test);
    out = read_meemum(dir=dir_test, filename=file_test);

    compare_meemum_results(out, out_expected)       # test if they are sufficiently similar

end

@testset "Vertex klb691 test" begin
    testname = "klb691"

    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", vertex=true) == true # run file
    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", pssect=true) == true # create plot

    @test  filesize("$(testname)/$(testname).blk") ≈  filesize("$(testname)/output/$(testname).blk") rtol=0.02      # number optim.
    @test  filesize("$(testname)/$(testname).ps") ≈  filesize("$(testname)/output/$(testname).ps") rtol=0.07        # size plots
  
    # number of performed optimizations  
    minG          = extract_value("$(testname)/$(testname).tim", "SQP G evaluations", " ")
    minG_expected = extract_value("$(testname)/output/$(testname).tim", "SQP G evaluations", " ")
    @test minG ≈ minG_expected rtol=1e-2

end


@testset "Vertex bl691 test" begin

    testname = "bl691"

    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", vertex=true) == true # run file
    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", pssect=true) == true # create plot

    @test  filesize("$(testname)/$(testname).blk") ≈  filesize("$(testname)/output/$(testname).blk") rtol=0.02      # number optim.
    @test  filesize("$(testname)/$(testname).ps") ≈  filesize("$(testname)/output/$(testname).ps") rtol=0.07        # size plots

    # number of performed optimizations  
    minG          = extract_value("$(testname)/$(testname).tim", "SQP G evaluations", " ")
    minG_expected = extract_value("$(testname)/output/$(testname).tim", "SQP G evaluations", " ")
    @test minG ≈ minG_expected rtol=1e-2    

end

@testset "Vertex weigang test" begin

    testname = "weigang"

    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", vertex=true) == true # run file
    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", pssect=true) == true # create plot

    @test  filesize("$(testname)/$(testname).blk") ≈  filesize("$(testname)/output/$(testname).blk") rtol=0.02      # number optim.
    @test  filesize("$(testname)/$(testname).ps") ≈  filesize("$(testname)/output/$(testname).ps") rtol=0.07        # size plots

    # number of performed optimizations  
    minG          = extract_value("$(testname)/$(testname).tim", "SQP G evaluations", " ")
    minG_expected = extract_value("$(testname)/output/$(testname).tim", "SQP G evaluations", " ")
    @test minG ≈ minG_expected rtol=1e-2    
end
