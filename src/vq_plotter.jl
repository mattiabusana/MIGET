using DataFrames, CairoMakie, Interpolations


function plot_defined_vq(qt, sdq1, sdq2, qmean1, qmean2, qratio, vdvt, qsqt; saving = false, name_saved_file = "defined_vq", file_extension = "png")


    qrest = qt * (1 - qsqt / 100)
    vqlo = 0.005
    vqhi = 100
    ncomp = 48
    slo = 0.001
    shi = 1000
    nsol = 11


    df_comps = DataFrame()
    df_comps.comparts = range(1, 49)
    b1 = 1
    b2 = 2
    b3 = 3
    b4 = 4
    b5 = 5
    b6 = 6
    b7 = 7
    b8 = 8
    b9 = 9
    b10 = 10
    b11 = 11

    delta = (log(vqhi / vqlo)) / (ncomp - 1)
    df_comps.vq_comparts = @. vqlo * (exp(delta * (df_comps.comparts - 1)))

    d11 = @. log(qmean1) .- log(df_comps.vq_comparts)
    d12 = @. log(qmean2) .- log(df_comps.vq_comparts)
    rexp1 = @. -0.5 * d11 * d11 / (sdq1 * sdq1)
    rexp2 = @. -0.5 * d12 * d12 / (sdq2 * sdq2)

    df_comps.main_mode_bf = exp.(rexp1)
    df_comps.secondary_mode_bf = qratio .* exp.(rexp2)
    df_comps.bf_combined = df_comps.main_mode_bf .+ df_comps.secondary_mode_bf
    sumq = sum(df_comps.bf_combined)

    df_comps.perfusion_normalized = qrest .* df_comps.bf_combined ./ sumq
    df_comps.ventilation_normalized = df_comps.perfusion_normalized .* df_comps.vq_comparts


    va = sum(df_comps.ventilation_normalized)
    ve = va / (1 - vdvt / 100)
    ########################


    df_gas = DataFrame()
    df_gas.gases = range(1, 11)

    deltas = (log(shi / slo)) / (nsol - 1)
    df_gas.partition_coefficients = @. slo * (exp(deltas * (df_gas.gases - 1)))

    df_gas.retentions_no_shunt = repeat([NaN], 11)
    df_gas.excretions_no_deadspace = repeat([NaN], 11)

    a1 = df_gas[b1, "partition_coefficients"]
    df_gas[b1, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a1 ./ (a1 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b1, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a1 ./ (a1 .+ df_comps.vq_comparts)) ./ va


    a2 = df_gas[b2, "partition_coefficients"]
    df_gas[b2, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a2 ./ (a2 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b2, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a2 ./ (a2 .+ df_comps.vq_comparts)) ./ va


    a3 = df_gas[b3, "partition_coefficients"]
    df_gas[b3, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a3 ./ (a3 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b3, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a3 ./ (a3 .+ df_comps.vq_comparts)) ./ va


    a4 = df_gas[b4, "partition_coefficients"]
    df_gas[b4, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a4 ./ (a4 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b4, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a4 ./ (a4 .+ df_comps.vq_comparts)) ./ va


    a4 = df_gas[b4, "partition_coefficients"]
    df_gas[b4, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a4 ./ (a4 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b4, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a4 ./ (a4 .+ df_comps.vq_comparts)) ./ va

    a4 = df_gas[b4, "partition_coefficients"]
    df_gas[b4, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a4 ./ (a4 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b4, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a4 ./ (a4 .+ df_comps.vq_comparts)) ./ va

    a5 = df_gas[b5, "partition_coefficients"]
    df_gas[b5, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a5 ./ (a5 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b5, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a5 ./ (a5 .+ df_comps.vq_comparts)) ./ va


    a6 = df_gas[b6, "partition_coefficients"]
    df_gas[b6, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a6 ./ (a6 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b6, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a6 ./ (a6 .+ df_comps.vq_comparts)) ./ va


    a7 = df_gas[b7, "partition_coefficients"]
    df_gas[b7, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a7 ./ (a7 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b7, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a7 ./ (a7 .+ df_comps.vq_comparts)) ./ va


    a8 = df_gas[b8, "partition_coefficients"]
    df_gas[b8, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a8 ./ (a8 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b8, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a8 ./ (a8 .+ df_comps.vq_comparts)) ./ va

    a9 = df_gas[b9, "partition_coefficients"]
    df_gas[b9, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a9 ./ (a9 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b9, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a9 ./ (a9 .+ df_comps.vq_comparts)) ./ va


    a10 = df_gas[b10, "partition_coefficients"]
    df_gas[b10, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a10 ./ (a10 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b10, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a10 ./ (a10 .+ df_comps.vq_comparts)) ./ va

    a10 = df_gas[b10, "partition_coefficients"]
    df_gas[b10, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a10 ./ (a10 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b10, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a10 ./ (a10 .+ df_comps.vq_comparts)) ./ va

    a11 = df_gas[b11, "partition_coefficients"]
    df_gas[b11, "retentions_no_shunt"] = sum(df_comps.perfusion_normalized .* a11 ./ (a11 .+ df_comps.vq_comparts)) ./ qrest
    df_gas[b11, "excretions_no_deadspace"] = sum(df_comps.ventilation_normalized .* a11 ./ (a11 .+ df_comps.vq_comparts)) ./ va


    df_gas.retentions_shunt = @. qsqt / 100 + (1 - qsqt / 100) * df_gas.retentions_no_shunt
    df_gas.excretions_deadspace = @. (1 - vdvt / 100) * df_gas.excretions_no_deadspace
    df_gas.retentions_homogeneous = @. df_gas.partition_coefficients / (df_gas.partition_coefficients + va / qt)
    df_gas.excretions_homogeneous = @. (1 - vdvt / 100) * df_gas.partition_coefficients / (df_gas.partition_coefficients + va / qt)

    df_shunt = DataFrame()
    df_shunt.vq_shunt = [0.0016, 0.0016]
    df_shunt.bf_shunt = [(0.01 * qsqt * qt), 0]
    df_shunt.ventilation_shunt = repeat([0], 2)

    df_shunt.vq_dead = [200, 200]
    df_shunt.bf_dead = [0.01 * vdvt * ve, 0]

    new_x = range(0.005, 100, length = 100000)
    itp_perf = interpolate(df_comps.vq_comparts, df_comps.perfusion_normalized, FritschCarlsonMonotonicInterpolation())
    itp_vent = interpolate(df_comps.vq_comparts, df_comps.ventilation_normalized, FritschCarlsonMonotonicInterpolation())
    CairoMakie.activate!(type = "png")
    #Makie.inline!(false)


    fig = Figure(resolution = (1600, 600))
    ax1 = Axis(fig[1, 1], xlabel = "Ventilation / Perfusion", ylabel = "Blood flow & ventilation (l/min)", xscale = log10, yminorticks = IntervalsBetween(3), yminorticksvisible = true, yminorgridvisible = true, xminorticksvisible = true,
        xminorticks = IntervalsBetween(8), xminorgridvisible = true, yticks = LinearTicks(4), xlabelsize = 18)
    scatter!(ax1, df_comps.vq_comparts, df_comps.perfusion_normalized, color = (:firebrick, 0.9), markersize = 10, label = nothing)
    scatter!(ax1, df_comps.vq_comparts, df_comps.ventilation_normalized, label = nothing, color = (:cornflowerblue, 0.9), markersize = 15, marker = '◼', markerstrokewidth = 0.1)
    lines!(ax1, new_x, itp_perf.(new_x), label = "Perfusion", color = :firebrick)
    lines!(ax1, new_x, itp_vent.(new_x), label = "Ventilation", color = :cornflowerblue, linestyle = :dash)

    if maximum(df_comps.perfusion_normalized) > maximum(df_comps.ventilation_normalized)
        whoismax = maximum(df_comps.perfusion_normalized)
    else
        whoismax = maximum(df_comps.ventilation_normalized)
    end


    if maximum(df_comps.perfusion_normalized) > df_shunt.bf_shunt[1]

        lines!(ax1, df_shunt.vq_shunt, df_shunt.bf_shunt, linewidth = 0.5, color = :firebrick, label = nothing)
        scatter!(ax1, (df_shunt.vq_shunt[1], df_shunt.bf_shunt[1]), color = (:firebrick, 0.9), label = nothing, markersize = 10)
        text!("Shunt = $qsqt%", position = (df_shunt.vq_shunt[1] - 0.0002, df_shunt.bf_shunt[1] + 0.05), color = :firebrick, textsize = 15)

    else
        arrows!([df_shunt.vq_shunt[1]], [whoismax - 0.1], [0], [0.1], color = (:firebrick, 0.9))
        scatter!(ax1, (df_shunt.vq_shunt[1], whoismax - 0.1), color = (:firebrick, 0.9), label = nothing, markersize = 10)
        text!("Shunt = $qsqt%", position = (df_shunt.vq_shunt[1] - 0.0002, whoismax - 0.1), color = :firebrick, textsize = 15)

    end


    if maximum(df_comps.ventilation_normalized) > df_shunt.bf_dead[1]
        lines!(ax1, df_shunt.vq_dead, df_shunt.bf_dead, linewidth = 0.5, color = :cornflowerblue, label = nothing)
        scatter!(ax1, (df_shunt.vq_dead[1], df_shunt.bf_dead[1]), color = (:cornflowerblue, 0.9), label = nothing, markersize = 15, marker = '◼')
        text!("Deadspace = $vdvt%", position = (df_shunt.vq_dead[1] - 175, df_shunt.bf_dead[1] + 0.05), color = :cornflowerblue, textsize = 15)

    else
        arrows!([df_shunt.vq_dead[1]], [whoismax - 0.1], [0], [0.1], color = (:cornflowerblue, 0.9))
        scatter!(ax1, (df_shunt.vq_dead[1], whoismax - 0.1), color = (:cornflowerblue, 0.9), label = nothing, markersize = 15, marker = '◼')
        text!("Deadspace = $vdvt%", position = (df_shunt.vq_dead[1] - 175, whoismax - 0.1), color = :cornflowerblue, textsize = 15)

    end



    ax2 = Axis(fig[1, 2], xlabel = "Gas-blood partition coefficient", ylabel = "Retention and excretion", xscale = log10, xminorgridvisible = true, xminorticksvisible = true, xminorticks = IntervalsBetween(8),
        yminorticks = IntervalsBetween(3), yminorticksvisible = true, yminorgridvisible = true)

    lines!(ax2, df_gas.partition_coefficients, df_gas.retentions_homogeneous, xaxis = :log, color = :steelblue, linestyle = :dash, linewidth = 1, label = "Retentions, homogeneous")
    scatter!(ax2, df_gas.partition_coefficients, df_gas.retentions_homogeneous, markersize = 10, color = :steelblue, label = nothing)

    lines!(ax2, df_gas.partition_coefficients, df_gas.retentions_shunt, c = :steelblue, linestyle = :solid, linewidth = 2, label = "Retentions, measured")
    scatter!(ax2, df_gas.partition_coefficients, df_gas.retentions_shunt, markersize = 10, color = :steelblue, label = nothing)

    lines!(ax2, df_gas.partition_coefficients, df_gas.excretions_homogeneous, xaxis = :log, color = :darkorange, linestyle = :dash, linewidth = 1, label = "Excretions, homogeneous")
    scatter!(ax2, df_gas.partition_coefficients, df_gas.excretions_homogeneous, markersize = 10, color = :darkorange, label = nothing)

    lines!(ax2, df_gas.partition_coefficients, df_gas.excretions_deadspace, c = :darkorange, linestyle = :solid, linewidth = 2, label = "Excretions, measured")
    scatter!(ax2, df_gas.partition_coefficients, df_gas.excretions_deadspace, markersize = 10, color = :darkorange, label = nothing)

    display(fig)


    if saving == true
        name_file_final = name_saved_file * "." * file_extension
        save(name_file_final, fig, px_per_unit = 5)

    end



end


function plot_miget_outputx(path_file::String; saving = false, name_saved_file = "miget_output", file_extension = "png")


    f = readlines(path_file)

    comparts = []
    perfusions = []
    ventilations = []


    for line in f
        arr_line = split(line)
        push!(comparts, arr_line[1])
        push!(perfusions, arr_line[3])
        push!(ventilations, arr_line[2])
    end

    shunt = parse(Float32, split(f[1])[4])
    deadspace = parse(Float32, split(f[end])[5])

    comparts = parse.(Float32, comparts)
    perfusions = parse.(Float32, perfusions)
    ventilations = parse.(Float32, ventilations)

    deleteat!(comparts, 1)
    deleteat!(perfusions, 1)
    deleteat!(ventilations, 1)

    deleteat!(comparts, length(comparts))
    deleteat!(perfusions, length(perfusions))
    deleteat!(ventilations, length(ventilations))

    qs = shunt / (shunt + sum(perfusions))
    qs_100 = round(qs * 100, digits = 1)

    vd_vt = deadspace / (sum(ventilations) + deadspace)
    vd_vt_100 = round(vd_vt * 100, digits = 1)

    new_x = range(minimum(comparts), maximum(comparts), length = 100000)
    itp_perf = interpolate(comparts, perfusions, FritschCarlsonMonotonicInterpolation())
    itp_vent = interpolate(comparts, ventilations, FritschCarlsonMonotonicInterpolation())


    new_x = collect(new_x)
    perf_interpolated = convert.(Float32, itp_perf.(new_x))
    vent_interpolated = convert.(Float32, itp_vent.(new_x))


    CairoMakie.activate!(type = "png")


    fig = Figure(resolution = (800, 600))
    ax1 = Axis(fig[1, 1], xlabel = "Ventilation / Perfusion", ylabel = "Blood flow & ventilation (l/min)", xscale = log10, yminorticks = IntervalsBetween(3), yminorticksvisible = true, yminorgridvisible = true, xminorticksvisible = true,
        xminorticks = IntervalsBetween(8), xminorgridvisible = true, yticks = LinearTicks(4), xlabelsize = 18)
    scatter!(ax1, comparts, perfusions, color = (:firebrick, 0.9), markersize = 10)
    scatter!(ax1, comparts, ventilations, color = (:cornflowerblue, 0.9), markersize = 15, marker = '◼', markerstrokewidth = 0.1)
    lines!(ax1, new_x, perf_interpolated, label = "Perfusion", color = :firebrick)
    lines!(ax1, new_x, vent_interpolated, label = "Ventilation", color = :cornflowerblue, linestyle = :dash)


    if maximum(perfusions) > maximum(ventilations)
        whoismax = maximum(perfusions)
    else
        whoismax = maximum(ventilations)

    end


    min_range = Float32(0.0016)
    shunt = Float32(shunt)
    deadspace = Float32(deadspace)
    max_range = Float32(1000)


    if whoismax <= 0.8
        r = 0.14
    elseif whoismax < 1.5 && whoismax >= 0.5
        r = 0.18
    elseif whoismax < 2.5 && whoismax >= 1.5
        r =  0.22
    end


    if maximum(perfusions) > shunt

        lines!(ax1, [min_range, min_range], [0, shunt], linewidth = 0.5, color = :firebrick)
        scatter!(ax1, (min_range, shunt), color = (:firebrick, 0.9), markersize = 10)
        text!("Shunt = $qs_100%", position = (0.0016 - 0.0006, shunt + 0.05), color = :firebrick, textsize = 15)


    else
        arrows!([min_range], [whoismax - 0.1], [0], [0.1], color = (:firebrick, 0.9))
        scatter!(ax1, (min_range, whoismax - 0.1), color = (:firebrick, 0.9), markersize = 10)
        text!("Shunt = $qs_100%", position = (0.0016 - 0.0006, whoismax - 0.15), color = :firebrick, textsize = 15)

    end


    if maximum(ventilations) > deadspace
        lines!(ax1, [max_range], [deadspace], linewidth = 0.5, color = :cornflowerblue,)
        scatter!(ax1, (max_range, deadspace), color = (:cornflowerblue, 0.9), markersize = 15, marker = '◼')
        text!("Deadspace = $qs_100%", position = (1000 - 940, deadspace + 0.05), color = :cornflowerblue, textsize = 15)

    else
        arrows!([max_range], [whoismax - 0.1], [0], [0.1], color = (:cornflowerblue, 0.9))
        scatter!(ax1, (max_range, whoismax - 0.1), color = (:cornflowerblue, 0.9), markersize = 15, marker = '◼')
        text!("Deadspace = $vd_vt_100%", position = (1000 - 940, whoismax - r), color = :cornflowerblue, textsize = 15)

    end

    axislegend(ax1, position = :rb)

    display(fig)



    if saving == true
        name_file_final = name_saved_file * "." * file_extension
        save(name_file_final, fig, px_per_unit = 5)

    end


end


function plot_papv_pepv(dataset::DataFrame, row_selected, sample; corrected = true)
    
    λs = Vector(dataset[row_selected, [:λ_sf6_def, :λ_ethane_def, :λ_cyclo_def, :λ_iso_def, :λ_dth_def, :λ_acetone_def]])
    ref_ret = Vector(dataset[row_selected, [:papv_sf6_ref, :papv_ethane_ref, :papv_cyclo_ref, :papv_iso_ref, :papv_dth_ref, :papv_acetone_ref]])
    ref_exc = Vector(dataset[row_selected, [:pepv_sf6_ref, :pepv_ethane_ref, :pepv_cyclo_ref, :pepv_iso_ref, :pepv_dth_ref, :pepv_acetone_ref]])

    if sample == 1 && corrected == false
        retentions = Vector(dataset[row_selected, ["papv1_sf6", "papv1_ethane", "papv1_cyclo", "papv1_iso", "papv1_dth", "papv1_acetone"]])
        excretions = Vector(dataset[row_selected, ["pepv1_sf6", "pepv1_ethane", "pepv1_cyclo", "pepv1_iso", "pepv1_dth", "pepv1_acetone"]])

    elseif sample == 2 && corrected == false
        retentions = Vector(dataset[row_selected, ["papv2_sf6", "papv2_ethane", "papv2_cyclo", "papv2_iso", "papv2_dth", "papv2_acetone"]])
        excretions = Vector(dataset[row_selected, ["pepv2_sf6", "pepv2_ethane", "pepv2_cyclo", "pepv2_iso", "pepv2_dth", "pepv2_acetone"]])

    elseif sample == 1 && corrected == true
        retentions = Vector(dataset[row_selected, ["papv1_sf6", "papv1_ethane", "papv1_cyclo", "papv1_iso", "papv1_dth", "papv1_acetone_corrected"]])
        excretions = Vector(dataset[row_selected, ["pepv1_sf6", "pepv1_ethane", "pepv1_cyclo", "pepv1_iso", "pepv1_dth", "pepv1_acetone_corrected"]])

    elseif sample == 2 && corrected == true
        retentions = Vector(dataset[row_selected, ["papv2_sf6", "papv2_ethane", "papv2_cyclo", "papv2_iso", "papv2_dth", "papv2_acetone_corrected"]])
        excretions = Vector(dataset[row_selected, ["pepv2_sf6", "pepv2_ethane", "pepv2_cyclo", "pepv2_iso", "pepv2_dth", "pepv2_acetone_corrected"]])

    end


    CairoMakie.activate!(type = "png")


    fig = Figure(resolution = (1200, 500))

    #Retentions

    ax1 = Axis(fig[1, 1], xscale = log10, title = "Retentions", xlabel = "Blood:gas partition coefficient", ylabel = "Retention",
        yminorticks = IntervalsBetween(3), yminorticksvisible = true, yminorgridvisible = true, xminorticksvisible = true, xminorticks = IntervalsBetween(8), xminorgridvisible = true)
    scatter!(ax1, λs, ref_ret, markersize = 10, marker = '◼')
    lines!(ax1, λs, ref_ret, linestyle = :dash, label = "Reference", linewidth = 1)
    scatter!(ax1, λs, retentions, markersize = 10)
    lines!(ax1, λs, retentions, label = "Patient")

    axislegend(ax1, position = :rb)

    # Excretions

    ax2 = Axis(fig[1, 2], xscale = log10, title = "Excretions", xlabel = "Blood:gas partition coefficient", ylabel = "Excretion",
        yminorticks = IntervalsBetween(3), yminorticksvisible = true, yminorgridvisible = true, xminorticksvisible = true, xminorticks = IntervalsBetween(8), xminorgridvisible = true)
    scatter!(ax2, λs, ref_exc, markersize = 10, marker = '◼')
    lines!(ax2, λs, ref_exc, label = "Reference", linestyle = :dash, linewidth = 1)
    scatter!(ax2, λs, excretions, markersize = 10)
    lines!(ax2, λs, excretions, label = "Patient")

    axislegend(ax2, position = :rb)


    display(fig)


    return nothing


end


# The next function takes a list of pairs name : best rss sample and a directory where gas or flows files are stored
# and outputs a list of only the files that need to be plotted. 

function select_best_rss_file(tup, directory)

    list_keep = []

    list_files = (sort ∘ readdir)(directory)
    sorted_tup = sort(tup)

    chunked = collect(Iterators.partition(list_files, 2))

    for (tup_element, patient) in zip(sorted_tup, chunked)

        cond = split(tup_element[2], "a")[2]
        
        if occursin(cond, patient[1])
            push!(list_keep, patient[1])
        elseif occursin(cond, patient[2])
            push!(list_keep, patient[2])
        end

    end


    return list_keep
end


