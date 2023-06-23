using GTOC12
using Test
using SPICE
furnish_all_kernels()

# @testset "GTOC1.jl" begin  
#     # define the test set
#     @testset "Type stability tests" begin
#         # test functions for type stability using @typewarn
#         @testset "Functions" for (name, method) in collect(methods(MyModule))
#             @testset "Function $name" begin
#                 @typewarn typeof(method(rand(eltype(method.parameters)))) 
#             end
#         end

#         # # test types for type stability
#         # @testset "Types" begin
#         #     @test isconcretetype(MyModule.MyType)
#         #     @test isbits(MyModule.MyType)
#         # end
#     end

# end

@testset "naive_flyby keplerian utilities" begin
    @testset "RV2COE Sun-centered" begin
        bodies = ["Earth", "MARS BARYCENTER", "VENUS BARYCENTER"]
        epochs = [86400, 86400 * (365 * 21 + 119)]
        CB = GTOC1.default_CB_str
        for bdy in bodies, et in epochs
            x⃗_test = spkgeo(bodn2c(bdy), et, GTOC1.default_ref_frame, bodn2c(CB))[1]
            ta = GTOC1.RV2COE(x⃗=x⃗_test, μ_CB_or_CB_name=CB)
            spcout = SPICE.oscltx(x⃗_test, et, bodvrd(CB, "GM")[1])
            tb = spcout[[10, 2, 3, 4, 5, 9]]
            @test ta ≈ tb
        end
    end
    @testset "RV2COE Planet-centered" begin
        bdy = "Moon"
        CB = "Earth"
        epochs = [86400, 86400 * (365 * 21 + 119)]
        for et in epochs
            x⃗_bdy = spkgeo(bodn2c(bdy), et, GTOC1.default_ref_frame, GTOC1.default_CB_idx)[1]
            x⃗_CB = spkgeo(bodn2c(CB), et, GTOC1.default_ref_frame, GTOC1.default_CB_idx)[1]
            x⃗_test = x⃗_bdy - x⃗_CB
            ta = GTOC1.RV2COE(x⃗=x⃗_test, μ_CB_or_CB_name=CB)
            spcout = SPICE.oscltx(x⃗_test, et, bodvrd(CB, "GM")[1])
            tb = spcout[[10, 2, 3, 4, 5, 9]]
            @test ta ≈ tb
        end
    end
    @testset "COE2RV Sun-centered" begin
        bodies = ["Earth", "MARS BARYCENTER", "VENUS BARYCENTER"]
        epochs = [86400, 86400 * (365 * 21 + 119)]
        CB = GTOC1.default_CB_str
        for bdy in bodies, et in epochs
            x⃗_test = spkgeo(bodn2c(bdy), et, GTOC1.default_ref_frame, bodn2c(CB))[1]
            elts = GTOC1.RV2COE(x⃗=x⃗_test, μ_CB_or_CB_name=CB)
            rp = elts[1] * (1 - elts[2])
            M0 = GTOC1.mean_anom(x⃗=x⃗_test, μ_CB_or_CB_name=CB)
            ta = SPICE.conics([rp, elts[2:5]..., M0, et, bodvrd(CB, "GM")[1]], et)
            tb = GTOC1.COE2RV(COE=elts, μ_CB_or_CB_name=CB)
            @test ta ≈ tb
        end
    end
    @testset "COE2RV Planet-centered" begin
        bdy = "Moon"
        CB = "Earth"
        epochs = [86400, 86400 * (365 * 21 + 119)]
        for et in epochs
            x⃗_bdy = spkgeo(bodn2c(bdy), et, GTOC1.default_ref_frame, GTOC1.default_CB_idx)[1]
            x⃗_CB = spkgeo(bodn2c(CB), et, GTOC1.default_ref_frame, GTOC1.default_CB_idx)[1]
            x⃗_test = x⃗_bdy - x⃗_CB
            elts = GTOC1.RV2COE(x⃗=x⃗_test, μ_CB_or_CB_name=CB)
            rp = elts[1] * (1 - elts[2])
            M0 = GTOC1.mean_anom(x⃗=x⃗_test, μ_CB_or_CB_name=CB)
            ta = SPICE.conics([rp, elts[2:5]..., M0, et, bodvrd(CB, "GM")[1]], et)
            tb = GTOC1.COE2RV(COE=elts, μ_CB_or_CB_name=CB)
            @test ta ≈ tb
        end
    end
# @testset "GTOC1.jl" begin  
#     # define the test set
#     @testset "Type stability tests" begin
#         # test functions for type stability using @typewarn
#         @testset "Functions" for (name, method) in collect(methods(MyModule))
#             @testset "Function $name" begin
#                 @typewarn typeof(method(rand(eltype(method.parameters)))) 
#             end
#         end
        
#         # # test types for type stability
#         # @testset "Types" begin
#         #     @test isconcretetype(MyModule.MyType)
#         #     @test isbits(MyModule.MyType)
#         # end
#     end

end
