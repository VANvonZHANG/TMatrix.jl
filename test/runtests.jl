using TMatrix
using Test
using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(TMatrix)
end

@testset "TMatrix.jl" begin
    # include("test_special_functions.jl")
    # include("test_vswf.jl")
    # include("test_mie.jl")
    @test true
end
