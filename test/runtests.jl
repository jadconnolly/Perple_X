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

    # Compare scalars
    Compare_Fields = [:H, :S, :Cp, :variance, :spec_Cp, :P, :T]
    @testset "test field $field" for field in Compare_Fields
        @test out[field] ≈ out_expected[field]
    end

    # Compare fields
    Compare_Fields = (:H2O, :MgO, :Al2O3, :K2O, :CaO, :TiO2, :FeO, :O2, :Na2O, :SiO2)
    @testset "test chem. pot. $field" for field in Compare_Fields
        @test out.chem_pot[field] ≈ out_expected.chem_pot[field] rtol=1e-4
    end

    # Compare phase composition
    Compare_Fields = (Symbol("wt%"), Symbol("vol%"), Symbol("mol%"), :mol, :H2O, :MgO, :Al2O3, :K2O, :CaO, :TiO2, :FeO, :O2, :Na2O, :SiO2)
    for i=1:length(out.phase_comp)
     @testset "test $(out.phase_comp[i].Name), $field" for field in Compare_Fields
            @test out.phase_comp[i][field] ≈ out_expected.phase_comp[i][field] rtol=2e-3
        end
    end

end

@testset "Vertex klb691 test" begin
    out = run_test_pipe(dir="klb691", inputfile="klb691_input.txt", vertex=true)  
    @test out == true
end

