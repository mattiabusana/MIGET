using DataFrames, PyCall, Conda, XLSX, CSV

function import_miget_dataset(path_to_file::String)

    f = XLSX.readxlsx(path_to_file)
    sheet = XLSX.sheetnames(f)[1]
    df = XLSX.eachtablerow(f[sheet]) |> DataFrame

    return df

end

function export_fortran_report(study_name, sets, pb, body_t, temp_bath, coeff_var, solubility_o2, n_gas, n_comparts, smooth, min_vq, max_vq, coeff_sf6,
    coeff_ethane, coeff_cyclo, coeff_iso, coeff_dth, coeff_acetone, apeak_sf6, apeak_ethane, apeak_cyclo, apeak_iso, apeak_dth, apeak_acetone,
    gpeak_sf6, gpeak_ethane, gpeak_cyclo, gpeak_iso, gpeak_dth, gpeak_acetone, ve, qt, t_ve, gas_volume_vial, blood_volume_vial, heparin_volume_vial, gas_volume_vial_mix,
    blood_volume_vial_mix, heparin_volume_vial_mix, hb, hct, vo2, vco2, tol, fio2, fico2, p50, po2, pco2, ph)


    Conda.add("fortranformat")


    lines = Vector{String}(undef, 13)
    ff = pyimport("fortranformat")

    lines = []

    # We save each single line in a string, we add each line to a list and we finally write line by line to a .RAW file

    l1_formatter = ff.FortranRecordWriter("(60A)")
    l1 = l1_formatter.write([study_name])

    l2_formatter = ff.FortranRecordWriter("(I4, F7.1, F6.1, F6.1, F8.2, F8.4)")
    l2 = l2_formatter.write([sets, pb, body_t, temp_bath, coeff_var, solubility_o2])



    l3_formatter = ff.FortranRecordWriter("(I4, I4, F6.1, F7.3, F8.2)")
    l3 = l3_formatter.write([n_gas, n_comparts, smooth, min_vq, max_vq])

    l4_formatter = ff.FortranRecordWriter("(6E12.3)")
    l4 = l4_formatter.write([coeff_sf6, coeff_ethane, coeff_cyclo, coeff_iso, coeff_dth, coeff_acetone])

    l4 = " " * l4

    l5_formatter = ff.FortranRecordWriter("(6F7.1)")
    l5 = l5_formatter.write([apeak_sf6, apeak_ethane, apeak_cyclo, apeak_iso, apeak_dth, apeak_acetone])

    l6 = "    1.     1.     1.     1.     1.     1."

    l7 = l5_formatter.write([gpeak_sf6, gpeak_ethane, gpeak_cyclo, gpeak_iso, gpeak_dth, gpeak_acetone])

    l8 = l6

    l9 = "    0.0    0.0    0.0    0.0    0.0    0.0"

    l10 = l6

    l11_formatter = ff.FortranRecordWriter("(F7.2, F6.2, F7.1, F6.1, F6.1, F7.2, F7.2, F7.2, F7.2, F7.2)")
    l11 = l11_formatter.write([ve, qt, pb, body_t, t_ve, gas_volume_vial, blood_volume_vial, heparin_volume_vial, gas_volume_vial_mix, blood_volume_vial_mix, heparin_volume_vial_mix])
    l11 = replace(l11, "\n" => "")

    l12_formatter = ff.FortranRecordWriter("(F6.1, F6.1, F7.1, F7.1, F9.2, F7.4,F6.4, F6.1)")
    l12 = l12_formatter.write([hb, hct, vo2, vco2, tol, fio2, fico2, p50])

    l13_formatter = ff.FortranRecordWriter("(F7.1, F7.1, F6.2)")
    l13 = l13_formatter.write([po2, pco2, ph])

    push!(lines, l1)
    push!(lines, l2)
    push!(lines, l3)
    push!(lines, l4)
    push!(lines, l5)
    push!(lines, l6)
    push!(lines, l7)
    push!(lines, l8)
    push!(lines, l9)
    push!(lines, l10)
    push!(lines, l11)
    push!(lines, l12)
    push!(lines, l13)

    return lines, study_name


end


function write_fortran_report(fortran_multistring, file_name)

    name_with_extension = file_name * ".RAW"

    open(name_with_extension, "w") do io

        for l in fortran_multistring

            new_line = l * "\n"
            write(io, new_line)

        end

    end

    return nothing

end


function write_input_short(dataset::DataFrame, record_n; default_pc = false, scaling_peak_factor = 1e-4)

    get_line = dataset[record_n, :]
    study_name = get_line.surname

    sets = 1.0
    pb = 752
    coeff_var = 0.03
    temp_bath = 37.5
    solubility_o2 = 0.003
    n_gas = 6
    n_comparts = 50
    smooth = 40
    min_vq = 0.005
    max_vq = 100
    tol = 99000
    fico2 = 0
    p50 = 28.3

    body_t = get_line.temp
    t_ve = body_t


    # Choose if to use dafault or measured partition coefficients


    if default_pc == true
        coeff_sf6 = get_line.λ_sf6_def
        coeff_ethane = get_line.λ_ethane_def
        coeff_cyclo = get_line.λ_cyclo_def
        coeff_iso = get_line.λ_iso_def
        coeff_dth = get_line.λ_dth_def
        coeff_acetone = get_line.λ_acetone_def
    else
        coeff_sf6 = get_line.λ_sf6_measured
        coeff_ethane = get_line.λ_ethane_measured
        coeff_cyclo = get_line.λ_cyclo_measured
        coeff_iso = get_line.λ_iso_measured
        coeff_dth = get_line.λ_dth_measured
        coeff_acetone = get_line.λ_acetone_measured
    end
    # The MIGET measurements are always in double, i. e., we always have two samples of arterial blood and two samples of expired gas
    # Here we first load the set 1
    # The scaling factor is to comply with the original Fortran format used in SHORT. 

    apeak_sf6 = get_line.a1_sf6_area * scaling_peak_factor
    apeak_ethane = get_line.a1_ethane_area * scaling_peak_factor
    apeak_cyclo = get_line.a1_cyclo_area * scaling_peak_factor
    apeak_iso = get_line.a1_iso_area * scaling_peak_factor
    apeak_dth = get_line.a1_dth_area * scaling_peak_factor
    apeak_acetone = get_line.a1_acetone_area * scaling_peak_factor

    gpeak_sf6 = get_line.g1_sf6_area * scaling_peak_factor
    gpeak_ethane = get_line.g1_ethane_area * scaling_peak_factor
    gpeak_cyclo = get_line.g1_cyclo_area * scaling_peak_factor
    gpeak_iso = get_line.g1_iso_area * scaling_peak_factor
    gpeak_dth = get_line.g1_dth_area * scaling_peak_factor
    gpeak_acetone = get_line.g1_acetone_area * scaling_peak_factor

    ve = get_line.ve_a1
    qt = get_line.co_recalculated_a1

    gas_volume_vial = get_line.vial_a1_gas_volume
    blood_volume_vial = get_line.vial_a1_blood_volume
    heparin_volume_vial = 0

    gas_volume_vial_mix = 0
    blood_volume_vial_mix = 0
    heparin_volume_vial_mix = 0

    hb = get_line.art_hb_coox
    hct = get_line.art_hct
    vo2 = get_line.vo2
    vco2 = get_line.vco2
    fio2 = get_line.fio2
    po2 = get_line.art_po2
    pco2 = get_line.art_pco2
    ph = get_line.art_ph

    file_1_name = study_name * "_a1"

    # We write a .RAW file that will be fed into the SHORT Fortran Program written by P. D. Wagner

    raw_fortran_lines_1, _ = export_fortran_report(file_1_name, sets, pb, body_t, temp_bath, coeff_var, solubility_o2, n_gas, n_comparts, smooth, min_vq, max_vq, coeff_sf6,
        coeff_ethane, coeff_cyclo, coeff_iso, coeff_dth, coeff_acetone, apeak_sf6, apeak_ethane, apeak_cyclo, apeak_iso, apeak_dth, apeak_acetone,
        gpeak_sf6, gpeak_ethane, gpeak_cyclo, gpeak_iso, gpeak_dth, gpeak_acetone, ve, qt, t_ve, gas_volume_vial, blood_volume_vial, heparin_volume_vial, gas_volume_vial_mix,
        blood_volume_vial_mix, heparin_volume_vial_mix, hb, hct, vo2, vco2, tol, fio2, fico2, p50, po2, pco2, ph)

    write_fortran_report(raw_fortran_lines_1, file_1_name)


    # We load also the set 2

    apeak_sf6 = get_line.a2_sf6_area * scaling_peak_factor
    apeak_ethane = get_line.a2_ethane_area * scaling_peak_factor
    apeak_cyclo = get_line.a2_cyclo_area * scaling_peak_factor
    apeak_iso = get_line.a2_iso_area * scaling_peak_factor
    apeak_dth = get_line.a2_dth_area * scaling_peak_factor
    apeak_acetone = get_line.a2_acetone_area * scaling_peak_factor

    gpeak_sf6 = get_line.g2_sf6_area * scaling_peak_factor
    gpeak_ethane = get_line.g2_ethane_area * scaling_peak_factor
    gpeak_cyclo = get_line.g2_cyclo_area * scaling_peak_factor
    gpeak_iso = get_line.g2_iso_area * scaling_peak_factor
    gpeak_dth = get_line.g2_dth_area * scaling_peak_factor
    gpeak_acetone = get_line.g2_acetone_area * scaling_peak_factor

    ve = get_line.ve_a2
    qt = get_line.co_recalculated_a2

    gas_volume_vial = get_line.vial_a2_gas_volume
    blood_volume_vial = get_line.vial_a2_blood_volume
    heparin_volume_vial = 0

    gas_volume_vial_mix = 0
    blood_volume_vial_mix = 0
    heparin_volume_vial_mix = 0

    file_2_name = study_name * "_a2"

    # We write the second set in a separate .RAW file 

    raw_fortran_lines_2, _ = export_fortran_report(file_2_name, sets, pb, body_t, temp_bath, coeff_var, solubility_o2, n_gas, n_comparts, smooth, min_vq, max_vq, coeff_sf6,
        coeff_ethane, coeff_cyclo, coeff_iso, coeff_dth, coeff_acetone, apeak_sf6, apeak_ethane, apeak_cyclo, apeak_iso, apeak_dth, apeak_acetone,
        gpeak_sf6, gpeak_ethane, gpeak_cyclo, gpeak_iso, gpeak_dth, gpeak_acetone, ve, qt, t_ve, gas_volume_vial, blood_volume_vial, heparin_volume_vial, gas_volume_vial_mix,
        blood_volume_vial_mix, heparin_volume_vial_mix, hb, hct, vo2, vco2, tol, fio2, fico2, p50, po2, pco2, ph)

    write_fortran_report(raw_fortran_lines_2, file_2_name)

    dir = pwd()
    list_files = readdir(dir)


    format1 = file_1_name * ".RAW"
    format2 = file_2_name * ".RAW"

    trues = [occursin(format1, x) || occursin(format2, x) for x in list_files]
    number_trues = sum(trues)

    if number_trues == 2
        println("------------------------")
        println("Export of .RAW files was successful!")
        println("------------------------")
    else
        println("------------------------")
        println("Something went wrong...")
        println("------------------------")

    end


    return nothing



end


function import_ct_scans(path_to_file::String)

    data_cts = CSV.read(path_to_file, DataFrame, decimal = ',')

    names_data = names(data_cts)
    new_names = replace.(names_data, "." => "_")
    rename!(data_cts, new_names)

    data_cts.id = lowercase.(data_cts.id)

    data_cts = data_cts[!, 1:end-1]

    return data_cts

end


function import_ct_scans!(dataset::DataFrame, path_to_file::String)

    data_cts = CSV.read(path_to_file, DataFrame, decimal = ',')
    names_data = names(data_cts)
    new_names = replace.(names_data, "." => "_")
    rename!(data_cts, new_names)

    data_cts.id = lowercase.(data_cts.id)

    data_cts = data_cts[!, 1:end-1]

    leftjoin!(dataset, data_cts, on = [:surname => :id])


    return nothing


end


function import_output_vqbohr(path_to_file::String)

        f = readlines(path_to_file)

        global counter_non_empty = []


        global title = ""
        global rss = 0
        global bf_shunt = 1.0
        global bf_01_001 = 1.0
        global vent_01_001 = 1.0
        global bf_001_01 = 1.0
        global vent_001_01 = 1.0
        global bf_01_1 = 1.0
        global vent_01_1 = 1.0
        global bf_1_10 = 1.0
        global vent_1_10 = 1.0
        global bf_10_100 = 1.0
        global vent_10_100 = 1.0
        global vent_dead = 1.0
        global qmean = 1.0
        global sdq = 1.0
        global skewq = 1.0
        global meanva = 1.0
        global sdva = 1.0
        global skewva = 1.0
        global re_star_sf6 = 1.0
        global re_star_ethane = 1.0
        global re_star_cyclo = 1.0
        global re_star_iso = 1.0
        global re_star_dth = 1.0
        global re_star_acetone = 1.0
        global predicted_measured_po2_diff = 1.0


    for (i, line) in enumerate(f)

        

            if i == 6


                title = split(line, ":")[end]
                title = strip(title)


            elseif i == 22

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end
            
            elseif i == 23

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end
            elseif i == 24

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end
            elseif i == 25

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end
            elseif i == 26

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end
            elseif i == 27

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end
            elseif i == 28

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end
            elseif i == 29

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end
            elseif i == 30

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end

            elseif i == 31

                check_empty = isempty(line)
                if check_empty 
                    push!(counter_non_empty, i)
                end

        elseif i == 32

            check_empty = isempty(line)
            if check_empty 
                push!(counter_non_empty, i)
            end

        

        elseif i == 33

            check_empty = isempty(line)
            if check_empty 
                push!(counter_non_empty, i)
            end


        elseif i == 34

            check_empty = isempty(line)
            if check_empty 
                push!(counter_non_empty, i)
            end


        elseif i == 35

            check_empty = isempty(line)
            if check_empty 
                push!(counter_non_empty, i)
            end

            end

    end


    for (i, line) in enumerate(f)

        if i == 35 + counter_non_empty[1] - 22

            rss = parse(Float64, split(line)[end])

        elseif i == 40 + counter_non_empty[1] - 22

            bf_shunt = parse(Float64, split(line)[4])

        elseif i == 41 + counter_non_empty[1] - 22


            bf_01_001 = parse(Float64, split(line)[6])
            vent_01_001 = parse(Float64, split(line)[7])


        elseif i == 42 + counter_non_empty[1] - 22


            bf_001_01 = parse(Float64, split(line)[6])
            vent_001_01 = parse(Float64, split(line)[7])

        elseif i == 43 + counter_non_empty[1] - 22

            bf_01_1 = parse(Float64, split(line)[6])
            vent_01_1 = parse(Float64, split(line)[7])

        elseif i == 44 + counter_non_empty[1] - 22

            bf_1_10 = parse(Float64, split(line)[6])
            vent_1_10 = parse(Float64, split(line)[7])

        elseif i == 45 + counter_non_empty[1] - 22

            bf_10_100 = parse(Float64, split(line)[6])
            vent_10_100 = parse(Float64, split(line)[7])

        elseif i == 46 + counter_non_empty[1] - 22

            vent_dead = parse(Float64, split(line)[5])

        elseif i == 49 +counter_non_empty[1] - 22

            qmean = parse(Float64, split(line)[7])

        elseif i == 50 + counter_non_empty[1] - 22

            sdq = parse(Float64, split(line)[8])

        elseif i == 51 + counter_non_empty[1] - 22

            skewq = parse(Float64, split(line)[8])

        elseif i == 53 + counter_non_empty[1] - 22

            meanva = parse(Float64, split(line)[6])
        elseif i == 54 + counter_non_empty[1] - 22

            sdva = parse(Float64, split(line)[7])
        elseif i == 55 + counter_non_empty[1] - 22
            skewva = parse(Float64, split(line)[7])


        elseif i == 59 + counter_non_empty[1] - 22
            re_star_sf6 = parse(Float64, split(line)[10])
        elseif i == 60 + counter_non_empty[1] - 22
            re_star_ethane = parse(Float64, split(line)[10])
        elseif i == 61 + counter_non_empty[1] - 22
            re_star_cyclo = parse(Float64, split(line)[10])
        elseif i == 62 + counter_non_empty[1] - 22
            re_star_iso = parse(Float64, split(line)[10])
        elseif i == 63 + counter_non_empty[1] - 22
            re_star_dth = parse(Float64, split(line)[10])
        elseif i == 64 + counter_non_empty[1] - 22
            re_star_acetone = parse(Float64, split(line)[10])


        elseif i == 279 + counter_non_empty[1] - 22
            predicted_measured_po2_diff = parse(Float64, split(line)[7])


        end


    end

        pre = Dict("title" => title,
            "rss" => rss,
            "bf_shunt" => bf_shunt,
            "bf_01_001" => bf_01_001,
            "vent_01_001" => vent_01_001,
            "bf_001_01" => bf_001_01,
            "vent_001_01" => vent_001_01,
            "bf_01_1" => bf_01_1,
            "vent_01_1" => vent_01_1,
            "bf_1_10" => bf_1_10,
            "vent_1_10" => vent_1_10,
            "bf_10_100" => bf_10_100,
            "vent_10_100" => vent_10_100,
            "vent_dead" => vent_dead,
            "qmean" => qmean,
            "sdq" => sdq,
            "skewq" => skewq,
            "meanva" => meanva,
            "sdva" => sdva,
            "skewva" => skewva,
            "re_star_sf6" => re_star_sf6,
            "re_star_ethane" => re_star_ethane,
            "re_star_cyclo" => re_star_cyclo,
            "re_star_iso" => re_star_iso,
            "re_star_dth" => re_star_dth,
            "re_star_acetone" => re_star_acetone,
            "predicted_measured_po2_diff" => predicted_measured_po2_diff)


        surname = split(title, "_")[1]

        delete!(pre, "title")

        function rename_dict(name_dict, suffix)

            new_name = name_dict * "_" * suffix
            return new_name

        end


        names_dict = keys(pre)
        vals = values(pre)


        if contains(title, "a1")
            new_names = rename_dict.(names_dict, "a1")
            new_dict = Dict(zip(new_names, vals))
            push!(new_dict, "sample" => "a1")
        else
            new_names = rename_dict.(names_dict, "a2")
            new_dict = Dict(zip(new_names, vals))
            push!(new_dict, "sample" => "a2")

        end



        d = DataFrame(new_dict)

        d.surname .= surname

    return d

end


function import_batch_vqbohr(path_to_dir::String)

    list_bohr = readdir(path_to_dir)


    #=

    files_to_keep = []

    for f in list_bohr
        if contains(f, "BOHR") 
            new_name = r * f
            push!(files_to_keep, new_name)

        end

    end

    =#

    extended_names = joinpath.(path_to_dir, list_bohr)

    df_to_join = import_output_vqbohr.(extended_names)



    function discriminate_a1a2(list_dfs)

        a1s = []
        a2s = []

        for d in list_dfs

            if d.sample[1] == "a1"
                push!(a1s, d)
            elseif d.sample[1] == "a2"
                push!(a2s, d)

            end
        end

        return a1s, a2s
    end


    a1s, a2s = discriminate_a1a2(df_to_join)


    final1 = a1s[1]
    final2 = a2s[1]

    if length(a2s) > 1

        for twos in a2s[2:end]
            final2 = vcat(final2, twos)
        end
    else
    end

    if length(a1s) > 1

        for ones in a1s[2:end]
            final1 = vcat(final1, ones)
        end
    else
    end

    select!(final1, Not(:sample))
    select!(final2, Not(:sample))

    final_merged = outerjoin(final1, final2, on = :surname)

    function mean_cols(val1, val2)

        if ismissing(val1) && !ismissing(val2)
            r = val2
        elseif !ismissing(val1) && ismissing(val2)
            r = val1
        elseif ismissing(val1) && ismissing(val2)
            r = missing
        elseif !ismissing(val1) && !ismissing(val2)
            r = (val1 + val2) / 2

        end
        return r

    end

    transform!(final_merged, [:rss_a1, :rss_a2] => ByRow(mean_cols) => :mean_rss)
    transform!(final_merged, [:bf_shunt_a1, :bf_shunt_a2] => ByRow(mean_cols) => :mean_bf_shunt)

    transform!(final_merged, [:bf_01_001_a1, :bf_01_001_a2] => ByRow(mean_cols) => :mean_bf_01_001)
    transform!(final_merged, [:vent_01_001_a1, :vent_01_001_a2] => ByRow(mean_cols) => :mean_vent_01_001)
    transform!(final_merged, [:bf_001_01_a1, :bf_001_01_a2] => ByRow(mean_cols) => :mean_bf_001_01)
    transform!(final_merged, [:vent_001_01_a1, :vent_001_01_a2] => ByRow(mean_cols) => :mean_vent_001_01)
    transform!(final_merged, [:bf_01_1_a1, :bf_01_1_a2] => ByRow(mean_cols) => :mean_bf_01_1_01)
    transform!(final_merged, [:vent_01_1_a1, :vent_01_1_a2] => ByRow(mean_cols) => :mean_vent_01_1)
    transform!(final_merged, [:bf_1_10_a1, :bf_1_10_a2] => ByRow(mean_cols) => :mean_bf_1_10)
    transform!(final_merged, [:vent_1_10_a1, :vent_1_10_a2] => ByRow(mean_cols) => :mean_vent_1_10)
    transform!(final_merged, [:bf_10_100_a1, :bf_10_100_a2] => ByRow(mean_cols) => :mean_bf_10_100)
    transform!(final_merged, [:vent_10_100_a1, :vent_10_100_a2] => ByRow(mean_cols) => :mean_vent_10_100)

    transform!(final_merged, [:vent_dead_a1, :vent_dead_a2] => ByRow(mean_cols) => :mean_vent_dead)
    transform!(final_merged, [:qmean_a1, :qmean_a2] => ByRow(mean_cols) => :mean_qmean)
    transform!(final_merged, [:sdq_a1, :sdq_a2] => ByRow(mean_cols) => :mean_sdq)

    transform!(final_merged, [:skewq_a1, :skewq_a2] => ByRow(mean_cols) => :mean_skewq)
    transform!(final_merged, [:meanva_a1, :meanva_a2] => ByRow(mean_cols) => :mean_meanva)
    transform!(final_merged, [:sdva_a1, :sdva_a2] => ByRow(mean_cols) => :mean_sdva)
    transform!(final_merged, [:skewva_a1, :skewva_a2] => ByRow(mean_cols) => :mean_skewva)
    transform!(final_merged, [:re_star_sf6_a1, :re_star_sf6_a2] => ByRow(mean_cols) => :mean_re_star_sf6)
    transform!(final_merged, [:re_star_ethane_a1, :re_star_ethane_a2] => ByRow(mean_cols) => :mean_re_star_ethane)
    transform!(final_merged, [:re_star_cyclo_a1, :re_star_cyclo_a2] => ByRow(mean_cols) => :mean_re_star_cyclo)
    transform!(final_merged, [:re_star_iso_a1, :re_star_iso_a2] => ByRow(mean_cols) => :mean_re_star_iso)
    transform!(final_merged, [:re_star_dth_a1, :re_star_dth_a2] => ByRow(mean_cols) => :mean_re_star_dth)
    transform!(final_merged, [:re_star_acetone_a1, :re_star_acetone_a2] => ByRow(mean_cols) => :mean_re_star_acetone)
    transform!(final_merged, [:predicted_measured_po2_diff_a1, :predicted_measured_po2_diff_a2] => ByRow(mean_cols) => :mean_predicted_measured_po2_diff)

    return final_merged


end


function import_batch_vqbohr!(dataset::DataFrame, path_to_dir::String)

    list_bohr = readdir(path_to_dir)


    #=

    files_to_keep = []

    for f in list_bohr
        if contains(f, "BOHR") 
            new_name = r * f
            push!(files_to_keep, new_name)

        end

    end

    =#

    extended_names = joinpath.(path_to_dir, list_bohr)

    df_to_join = import_output_vqbohr.(extended_names)



    function discriminate_a1a2(list_dfs)

        a1s = []
        a2s = []

        for d in list_dfs

            if d.sample[1] == "a1"
                push!(a1s, d)
            elseif d.sample[1] == "a2"
                push!(a2s, d)

            end
        end

        return a1s, a2s
    end


    a1s, a2s = discriminate_a1a2(df_to_join)


    final1 = a1s[1]
    final2 = a2s[1]

    if length(a2s) > 1

        for twos in a2s[2:end]
            final2 = vcat(final2, twos)
        end
    else
    end

    if length(a1s) > 1

        for ones in a1s[2:end]
            final1 = vcat(final1, ones)
        end
    else
    end

    select!(final1, Not(:sample))
    select!(final2, Not(:sample))

    final_merged = outerjoin(final1, final2, on = :surname)

    function mean_cols(val1, val2)

        if ismissing(val1) && !ismissing(val2)
            r = val2
        elseif !ismissing(val1) && ismissing(val2)
            r = val1
        elseif ismissing(val1) && ismissing(val2)
            r = missing
        elseif !ismissing(val1) && !ismissing(val2)
            r = (val1 + val2) / 2

        end
        return r

    end

    transform!(final_merged, [:rss_a1, :rss_a2] => ByRow(mean_cols) => :mean_rss)
    transform!(final_merged, [:bf_shunt_a1, :bf_shunt_a2] => ByRow(mean_cols) => :mean_bf_shunt)

    transform!(final_merged, [:bf_01_001_a1, :bf_01_001_a2] => ByRow(mean_cols) => :mean_bf_01_001)
    transform!(final_merged, [:vent_01_001_a1, :vent_01_001_a2] => ByRow(mean_cols) => :mean_vent_01_001)
    transform!(final_merged, [:bf_001_01_a1, :bf_001_01_a2] => ByRow(mean_cols) => :mean_bf_001_01)
    transform!(final_merged, [:vent_001_01_a1, :vent_001_01_a2] => ByRow(mean_cols) => :mean_vent_001_01)
    transform!(final_merged, [:bf_01_1_a1, :bf_01_1_a2] => ByRow(mean_cols) => :mean_bf_01_1_01)
    transform!(final_merged, [:vent_01_1_a1, :vent_01_1_a2] => ByRow(mean_cols) => :mean_vent_01_1)
    transform!(final_merged, [:bf_1_10_a1, :bf_1_10_a2] => ByRow(mean_cols) => :mean_bf_1_10)
    transform!(final_merged, [:vent_1_10_a1, :vent_1_10_a2] => ByRow(mean_cols) => :mean_vent_1_10)
    transform!(final_merged, [:bf_10_100_a1, :bf_10_100_a2] => ByRow(mean_cols) => :mean_bf_10_100)
    transform!(final_merged, [:vent_10_100_a1, :vent_10_100_a2] => ByRow(mean_cols) => :mean_vent_10_100)

    transform!(final_merged, [:vent_dead_a1, :vent_dead_a2] => ByRow(mean_cols) => :mean_vent_dead)
    transform!(final_merged, [:qmean_a1, :qmean_a2] => ByRow(mean_cols) => :mean_qmean)
    transform!(final_merged, [:sdq_a1, :sdq_a2] => ByRow(mean_cols) => :mean_sdq)

    transform!(final_merged, [:skewq_a1, :skewq_a2] => ByRow(mean_cols) => :mean_skewq)
    transform!(final_merged, [:meanva_a1, :meanva_a2] => ByRow(mean_cols) => :mean_meanva)
    transform!(final_merged, [:sdva_a1, :sdva_a2] => ByRow(mean_cols) => :mean_sdva)
    transform!(final_merged, [:skewva_a1, :skewva_a2] => ByRow(mean_cols) => :mean_skewva)
    transform!(final_merged, [:re_star_sf6_a1, :re_star_sf6_a2] => ByRow(mean_cols) => :mean_re_star_sf6)
    transform!(final_merged, [:re_star_ethane_a1, :re_star_ethane_a2] => ByRow(mean_cols) => :mean_re_star_ethane)
    transform!(final_merged, [:re_star_cyclo_a1, :re_star_cyclo_a2] => ByRow(mean_cols) => :mean_re_star_cyclo)
    transform!(final_merged, [:re_star_iso_a1, :re_star_iso_a2] => ByRow(mean_cols) => :mean_re_star_iso)
    transform!(final_merged, [:re_star_dth_a1, :re_star_dth_a2] => ByRow(mean_cols) => :mean_re_star_dth)
    transform!(final_merged, [:re_star_acetone_a1, :re_star_acetone_a2] => ByRow(mean_cols) => :mean_re_star_acetone)
    transform!(final_merged, [:predicted_measured_po2_diff_a1, :predicted_measured_po2_diff_a2] => ByRow(mean_cols) => :mean_predicted_measured_po2_diff)

    leftjoin!(dataset, final_merged, on = :surname)



    return nothing


end