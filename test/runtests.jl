using Test, JSON, ModeMatching, Match, DelimitedFiles, Plots, Base.Filesystem, Base.Threads

@testset "Comparison to comsol" begin
    tests = JSON.parsefile("comsol_tests.json")
    for (name, t) in tests
        println("Running $name...")
        n_TE = 20
        n_TM = 20
        max_m = 4
        f = t["freq"]
        initial_mode = t["mode_in"]
        comsol_file = t["comsol_file"]
        waveguides = []
        for w in t["setup"]
            if w["type"] == "RectangularWaveguide"
                push!(waveguides, RectangularWaveguide(w["x"], w["y"], w["z"], w["a"], w["b"],
                                                       w["length"], w["μ"], w["ε"], f))
            end
        end
        if initial_mode == "TE10"
            a_1 = vcat(1 / sqrt(power(waveguides[1], TEMode(1, 0))), zeros(n_TE + n_TM - 1))
        else
            throw("Unimplemented")
        end
        setup = WaveguideSetup(waveguides, n_TE, n_TM, max_m, a_1)
        (data_raw, header_raw) = readdlm(t["comsol_file"], ' ', AbstractString;
                                 header=true, comments=true, comment_char='%')
        coordinates = parse.(Float64, @view data_raw[:, 1:3])
        data = parse.(ComplexF64, @view data_raw[:, 4:end])
        header = @view header_raw[4:end]

        (t_ab, r_ab) = t_r_ab(setup)
        a_i = calc_a_i(setup, t_ab, r_ab)
        b_i = calc_b_i(setup, t_ab, r_ab);

        mses = Array{AbstractFloat}(undef, size(data))
        #max = 0
        sim_data = Array{Complex}(undef, size(data))
        @threads for (ri, r)=collect(enumerate(eachrow(coordinates)))
            (x, y, z) = r
            e = E(setup, a_i, b_i, x, y, z)
            h = H(setup, a_i, b_i, x, y, z)
            for (ci, expr) in enumerate(header)
                sim = @match expr begin
                    "ewfd.Ex" => e[1]
                    "ewfd.Ey" => e[2]
                    "ewfd.Ez" => e[3]
                    "ewfd.Hx" => h[1]
                    "ewfd.Hy" => h[2]
                    "ewfd.Hz" => h[3]
                end
                if isnan(sim)
                    sim = 0
                end
                if isnan(data[ri, ci])
                    data[ri, ci] = 0
                end
                #=if abs(data[ri, ci]) > max
                    max = abs(data[ri, ci])
                end=#
                # Conjugate is necessary because comsol uses exp(-iβz) for forward propagating waves, we use exp(iβz)
                mses[ri, ci] = abs2(sim - conj(data[ri, ci]))
                sim_data[ri, ci] = sim
            end
        end
        ind_yslice = 0.35e-2 .< coordinates[:, 2] .< 0.45e-2
        gr()
        xs = @view coordinates[ind_yslice, 1]
        zs = @view coordinates[ind_yslice, 3]
        mkpath("results/$name")
        for (ci, expr) in enumerate(header)
            es_com = @view data[ind_yslice, ci]
            es_sim = @view sim_data[ind_yslice, ci]
            lims = (minimum(abs.(es_com)), maximum(abs.(es_com)))
            p = plot(layout=grid(1, 2, widths=[0.40, 0.60]), size=(600, 1000),
                     shape=:rect, clim=lims, leg=true, cb=false)
            scatter!(p, xs, zs, mz=abs.(es_com), sp=1, ms=3, msw=0, msa=0,
                     shape=:rect, clim=lims, label="Comsol")
            scatter!(p, xs, zs, mz=abs.(es_sim), sp=2, ms=3, msw=0, msa=0,
                     shape=:rect, clim=lims, cb=true, label="Mode matching")
            #scatter!(p, [NaN], [NaN], mz=[NaN], framestyle=:none, grid=false, leg=:false,
            #         clim=lims, cb=true, sp=3)
            savefig(joinpath(["results", name, "$expr.svg"]))
        end
        println("MSEs:")
        for (ci, expr) in enumerate(header)
            println("$expr: $(sum(mses[:, ci]) / length(mses[:, ci]))")
        end
    end
end
