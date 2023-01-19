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
    out = run_test_pipe(dir="klb691", inputfile="klb691_input.txt", vertex=true)  # run file
    @test out == true

    out = run_test_pipe(dir="klb691", inputfile="klb691_input.txt", pssect=true)  # create plot
    @test out == true

    # check number of optimizations:
    @test  filesize("klb691/klb691.blk") ≈  filesize("klb691/output/klb691.blk") rtol=0.02  

    # check size of generated plot
    @test  filesize("klb691/klb691.ps") ≈  filesize("klb691/output/klb691.ps") rtol=0.07    

end


@testset "Vertex bl691 test" begin

    testname = "bl691"

    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", vertex=true) == true # run file
    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", pssect=true) == true # create plot

    @test  filesize("$(testname)/$(testname).blk") ≈  filesize("$(testname)/output/$(testname).blk") rtol=0.02      # number optim.
    @test  filesize("$(testname)/$(testname).ps") ≈  filesize("$(testname)/output/$(testname).ps") rtol=0.07        # size plots

end

@testset "Vertex weigang test" begin

    testname = "weigang"

    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", vertex=true) == true # run file
    @test run_test_pipe(dir=testname, inputfile="$(testname)_input.txt", pssect=true) == true # create plot

    @test  filesize("$(testname)/$(testname).blk") ≈  filesize("$(testname)/output/$(testname).blk") rtol=0.02      # number optim.
    @test  filesize("$(testname)/$(testname).ps") ≈  filesize("$(testname)/output/$(testname).ps") rtol=0.07        # size plots

end
