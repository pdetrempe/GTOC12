using GTOC1
using Test

@testset "GTOC1.jl" begin  
    # define the test set
    @testset "Type stability tests" begin
        # test functions for type stability using @typewarn
        @testset "Functions" for (name, method) in collect(methods(MyModule))
            @testset "Function $name" begin
                @typewarn typeof(method(rand(eltype(method.parameters)))) 
            end
        end
        
        # # test types for type stability
        # @testset "Types" begin
        #     @test isconcretetype(MyModule.MyType)
        #     @test isbits(MyModule.MyType)
        # end
    end

end
