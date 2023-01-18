using Test

include("utils.jl")

@testset "MEEMUM tests" begin
    @test run_meemum(dir="bl691_MEEMUM", inputfile="bl691_MEEMUM_input.txt", dir_meemum="../sources/meemum") == true
end
