@testset "Gauss–Lobatto" begin
    # Check error
    @test_throws DomainError gausslobatto(-1)
    @test_throws DomainError gausslobatto(0)
    @test_throws DomainError gausslobatto(1)

    n = 2
    x,w = gausslobatto(n)
    @test x[1] ≈ -1.
    @test x[2] ≈ 1.
    @test w[1] ≈ 1.
    @test w[2] ≈ 1.

    n = 3
    x,w = gausslobatto(n)
    @test x[1] ≈ -1.
    @test abs(x[2])<1e-15
    @test x[3] ≈ 1.
    @test w[1] ≈ 1/3
    @test w[2] ≈ 4/3
    @test w[3] ≈ 1/3

    tol = 1e-14
    n = 42
    x,w = gausslobatto(n)
    @test ( length(x) == n) && ( length(w) == n )
    @test abs(x[37] - 0.922259214258616) < tol
    @test abs(w[37] - 0.029306411216166) < tol
    @test ( x[1] == -1 && x[n] == 1 )

    @test dot( w,(x.^2)) ≈ 2/3
    @test dot( w,exp.(x)) ≈ exp(1)-exp(-1)
end