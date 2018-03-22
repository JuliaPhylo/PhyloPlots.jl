@testset "test of substitution model plot, 2 binary traits" begin

m3 = TwoBinaryTraitSubstitutionModel([2.0,1.2,1.1,2.2,1.0,3.1,2.0,1.1],
     ["carnivory", "noncarnivory", "wet", "dry"]);
@test_nowarn plot(m3)

end
