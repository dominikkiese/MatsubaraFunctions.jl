function get_timers() :: Nothing 

    to = TimerOutput()

    # interpolation on linear grid
    @timeit to "Linear interpolation" begin
        fg = MatsubaraGrid(1.0, 64, Fermion); nf = length(fg)
        bg = MatsubaraGrid(1.0, 64, Boson);   nb = length(bg)

        f1D = MatsubaraFunction((fg,), (10, 10), rand(nf, 10, 10))
        f2D = MatsubaraFunction((bg, fg), (10, 10), rand(nb, nf, 10, 10))
        f3D = MatsubaraFunction((bg, fg, fg), (10, 10), rand(nb, nf, nf, 10, 10))

        x   = rand(1 : 10), rand(1 : 10)
        w1D = fg[rand(1 : nf)]
        w2D = bg[rand(1 : nb)], fg[rand(1 : nf)]
        w3D = bg[rand(1 : nb)], fg[rand(1 : nf)], fg[rand(1 : nf)]

        for n in 1 : 100
            @timeit to "-> 1D" f1D(w1D, x...) 
            @timeit to "-> 2D" f2D(w2D, x...) 
            @timeit to "-> 3D" f3D(w3D, x...) 
        end 
    end

    # interpolation on coarse grid
    @timeit to "Coarse interpolation" begin
        fg = MatsubaraGrid(1.0, 64, 2.0, Fermion); nf = length(fg)
        bg = MatsubaraGrid(1.0, 64, 2.0, Boson);   nb = length(bg)

        f1D = MatsubaraFunction((fg,), (10, 10), rand(nf, 10, 10))
        f2D = MatsubaraFunction((bg, fg), (10, 10), rand(nb, nf, 10, 10))
        f3D = MatsubaraFunction((bg, fg, fg), (10, 10), rand(nb, nf, nf, 10, 10))

        x   = rand(1 : 10), rand(1 : 10)
        w1D = fg[rand(1 : nf)]
        w2D = bg[rand(1 : nb)], fg[rand(1 : nf)]
        w3D = bg[rand(1 : nb)], fg[rand(1 : nf)], fg[rand(1 : nf)]

        for n in 1 : 100
            @timeit to "-> 1D" f1D(w1D, x...) 
            @timeit to "-> 2D" f2D(w2D, x...) 
            @timeit to "-> 3D" f3D(w3D, x...) 
        end 
    end

    # extrapolation using quadratic model
    @timeit to "Extrapolation" begin
        fg = MatsubaraGrid(1.0, 512, Fermion)
        w  = 2.0 * fg[end]
        f  = MatsubaraFunction((fg,), (1,))

        for v in 1 : length(fg)
            f[v, 1] = 1.0 / (im * fg[v])
        end 

        for n in 1 : 100
            @timeit to "-> tail moments"  tail_moments(f, 1)
            @timeit to "-> extrapolation" f(w, 1; extrp = true)
        end 
    end

    # matsubara summation
    @timeit to "Summation" begin
        fg = MatsubaraGrid(1.0, 512, Fermion)
        f1 = MatsubaraFunction((fg,), (1,), Float64)
        f2 = MatsubaraFunction((fg,), (1,))

        for v in 1 : length(fg)
            f1[v, 1] = 1.0 / fg[v]
            f2[v, 1] = 1.0 / (im * fg[v])
        end 

        for n in 1 : 100
            @timeit to "-> Q = Float64"    sum(f1, 1)
            @timeit to "-> Q = ComplexF64" sum(f2, 1)
        end 
    end

    show(to)
    return nothing 
end