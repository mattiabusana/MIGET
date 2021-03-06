using DataFramesMeta, Statistics

function miget_physiology_calculations!(dataset :: DataFrame)


    function calculate_bsa(sex, weight, height)

        if sex == "f"   # DuBois Formula
            bsa = 0.007184 * (weight^0.425) * (height^0.725)
            return bsa
        else
            bsa = 0.007184 * (weight^0.425) * (height^0.725)
            return bsa
        end

    end


    function o2_content(hb, po2, saturation)

        content = 1.39 * hb * (saturation / 100) + 0.00304 * po2

    end



    # Cardiac output during MIGET measurements

    dataset.bsa_recalculated = @with(dataset, calculate_bsa.(:sex, :weight, :height))

    dataset.co_recalculated = @with(dataset, :ci .* :bsa_recalculated)
    dataset.co_recalculated_a1 = @with(dataset, :ci_a1 .* :bsa_recalculated)
    dataset.co_recalculated_a2 = @with(dataset, :ci_a2 .* :bsa_recalculated)

    # Physiological measurements

    dataset.cao2 = @with(dataset, o2_content.(:art_hb_coox, :art_po2, :art_so2_calc))
    dataset.cvo2 = @with(dataset, o2_content.(:art_hb_coox, :ven_po2, :ven_so2_calc))

    dataset.vo2 = @. @with(dataset, 10 * :co_recalculated_a1 * (:cao2 - :cvo2) * 1.22)  # Check it!
    dataset.r = @with(dataset, :vco2 ./ :vo2)
    dataset.palvo2 = @. @with(dataset, 713 * :fio2 - ((:art_pco2 * (1 - :fio2 * (1 - :r))) / :r))
    dataset.aado2 = @with(dataset, :palvo2 - :art_po2)

    dataset.cco2 = @with(dataset, o2_content.(:art_hb_coox, :palvo2, 100))
    dataset.qsqt = @. @with(dataset, ((:cco2 - :cao2) / (:cco2 - :cvo2)))
    dataset.va = @. @with(dataset, (863 * :vco2) / :art_pco2)
    dataset.vdvt = @. @with(dataset, (:ve * 1000 - :va) / (:ve * 1000))

    function calculate_ve_measurement(ve, rr, tidal)
        ###### TO BE ADDED ####

    end


    return nothing


end


function miget_gases_calculations!(dataset :: DataFrame)



    # Parameters for blood density and quantities in the vials

    dataset.heparin_density .= 1
    dataset.vial_volume .= 12

    dataset.pooled_blood_hep_syr = @. @with(dataset, :c1_blood + :a1_blood + :a2_blood)
    dataset.pooled_hep_syr = @. @with(dataset, :c1_heparin + :a1_heparin + :a2_heparin)
    dataset.pooled_tare = @. @with(dataset, :c1_dry + :a1_dry + :a2_dry)

    dataset.pooled_heparin_weight = @. @with(dataset, :pooled_hep_syr - :pooled_tare)
    dataset.pooled_heparin_volume = @. @with(dataset, :pooled_heparin_weight / :heparin_density)
    dataset.pooled_blood_weight = @. @with(dataset, :pooled_blood_hep_syr - :pooled_heparin_weight - :pooled_tare)
    dataset.pooled_blood_volume = @. @with(dataset, (:c1_tot_vol + :a1_tot_vol + :a2_tot_vol) - :pooled_heparin_volume)
    dataset.pooled_blood_density = @. @with(dataset, :pooled_blood_weight / :pooled_blood_volume)

    m = median(skipmissing(dataset.pooled_blood_density))

    transform!(dataset, [:pooled_blood_density] => ByRow(x -> ismissing(x) ? m :
                                                    x) => :pooled_blood_density)

    dataset.vial_a1_blood_volume = @. @with(dataset, (:vial_a1_blood - :vial_a1_dry) / :pooled_blood_density)
    dataset.vial_a2_blood_volume = @. @with(dataset, (:vial_a2_blood - :vial_a2_dry) / :pooled_blood_density)


    dataset.vial_a1_gas_volume = @. @with(dataset, :vial_volume - :vial_a1_blood_volume)
    dataset.vial_a2_gas_volume = @. @with(dataset, :vial_volume - :vial_a2_blood_volume)

    # Default partition coefficients

    dataset.??_acetone_def .= 317
    dataset.??_dth_def .= 12.3
    dataset.??_iso_def .= 1.45
    dataset.??_cyclo_def .= 0.749
    dataset.??_ethane_def .= 0.132
    dataset.??_sf6_def .= 0.00658


    ### ATTENZIONE AL COEFF DI PART: default or not?  ######


    # Equilibration of gases in the vial

    volume_correction(area, blood_volume, gas_volume, ??) = @. area * (blood_volume * ??+ gas_volume) / (blood_volume * ??)


    transform!(dataset, [:a1_sf6_area, :vial_a1_blood_volume, :vial_a1_gas_volume, :??_sf6_def] => volume_correction => :a1_sf6_area_volumecorrected)
    transform!(dataset, [:a2_sf6_area, :vial_a2_blood_volume, :vial_a1_gas_volume, :??_sf6_def] => volume_correction => :a2_sf6_area_volumecorrected)

    transform!(dataset, [:a1_ethane_area, :vial_a1_blood_volume, :vial_a1_gas_volume, :??_ethane_def] => volume_correction => :a1_ethane_area_volumecorrected)
    transform!(dataset, [:a2_ethane_area, :vial_a2_blood_volume, :vial_a1_gas_volume, :??_ethane_def] => volume_correction => :a2_ethane_area_volumecorrected)

    transform!(dataset, [:a1_cyclo_area, :vial_a1_blood_volume,:vial_a1_gas_volume, :??_cyclo_def] => volume_correction => :a1_cyclo_area_volumecorrected)
    transform!(dataset, [:a2_cyclo_area, :vial_a2_blood_volume, :vial_a1_gas_volume,:??_cyclo_def] => volume_correction => :a2_cyclo_area_volumecorrected)

    transform!(dataset, [:a1_iso_area, :vial_a1_blood_volume,:vial_a1_gas_volume, :??_iso_def] => volume_correction => :a1_iso_area_volumecorrected)
    transform!(dataset, [:a2_iso_area, :vial_a2_blood_volume,:vial_a1_gas_volume, :??_iso_def] => volume_correction => :a2_iso_area_volumecorrected)

    transform!(dataset, [:a1_dth_area, :vial_a1_blood_volume,:vial_a1_gas_volume, :??_dth_def] => volume_correction => :a1_dth_area_volumecorrected)
    transform!(dataset, [:a2_dth_area, :vial_a2_blood_volume,:vial_a1_gas_volume, :??_dth_def] => volume_correction => :a2_dth_area_volumecorrected)

    transform!(dataset, [:a1_acetone_area, :vial_a1_blood_volume,:vial_a1_gas_volume, :??_acetone_def] => volume_correction => :a1_acetone_area_volumecorrected)
    transform!(dataset, [:a2_acetone_area, :vial_a2_blood_volume,:vial_a1_gas_volume, :??_acetone_def] => volume_correction => :a2_acetone_area_volumecorrected)



    # Acetone loss correction


    dataset.g1_acetone_area_corrected = @. dataset.a1_acetone_area_volumecorrected * dataset.g1_dth_area / dataset.a1_dth_area_volumecorrected
    dataset.g2_acetone_area_corrected = @. dataset.a2_acetone_area_volumecorrected * dataset.g2_dth_area / dataset.a2_dth_area_volumecorrected


    # Pv calculation through the Fick principle

    pv_calculation(a_area, g_area, ve, co, ??) = @. a_area + ((g_area * ve) / (co * ??))

    transform!(dataset, [:a1_sf6_area_volumecorrected, :g1_sf6_area, :ve_a1, :co_recalculated_a1, :??_sf6_def] => pv_calculation => :v1_sf6_area_fick)
    transform!(dataset, [:a2_sf6_area_volumecorrected, :g2_sf6_area, :ve_a2, :co_recalculated_a2, :??_sf6_def] => pv_calculation => :v2_sf6_area_fick)

    transform!(dataset, [:a1_ethane_area_volumecorrected, :g1_ethane_area, :ve_a1, :co_recalculated_a1, :??_ethane_def] => pv_calculation => :v1_ethane_area_fick)
    transform!(dataset, [:a2_ethane_area_volumecorrected, :g2_ethane_area, :ve_a2, :co_recalculated_a2, :??_ethane_def] => pv_calculation => :v2_ethane_area_fick)

    transform!(dataset, [:a1_cyclo_area_volumecorrected, :g1_cyclo_area, :ve_a1, :co_recalculated_a1, :??_cyclo_def] => pv_calculation => :v1_cyclo_area_fick)
    transform!(dataset, [:a2_cyclo_area_volumecorrected, :g2_cyclo_area, :ve_a2, :co_recalculated_a2, :??_cyclo_def] => pv_calculation => :v2_cyclo_area_fick)

    transform!(dataset, [:a1_iso_area_volumecorrected, :g1_iso_area, :ve_a1, :co_recalculated_a1, :??_iso_def] => pv_calculation => :v1_iso_area_fick)
    transform!(dataset, [:a2_iso_area_volumecorrected, :g2_iso_area, :ve_a2, :co_recalculated_a2, :??_iso_def] => pv_calculation => :v2_iso_area_fick)

    transform!(dataset, [:a1_dth_area_volumecorrected, :g1_dth_area, :ve_a1, :co_recalculated_a1, :??_dth_def] => pv_calculation => :v1_dth_area_fick)
    transform!(dataset, [:a2_dth_area_volumecorrected, :g2_dth_area, :ve_a2, :co_recalculated_a2, :??_dth_def] => pv_calculation => :v2_dth_area_fick)

    transform!(dataset, [:a1_acetone_area_volumecorrected, :g1_acetone_area, :ve_a1, :co_recalculated_a1, :??_acetone_def] => pv_calculation => :v1_acetone_area_fick)
    transform!(dataset, [:a2_acetone_area_volumecorrected, :g2_acetone_area, :ve_a2, :co_recalculated_a2, :??_acetone_def] => pv_calculation => :v2_acetone_area_fick)

    transform!(dataset, [:a1_acetone_area_volumecorrected, :g1_acetone_area_corrected, :ve_a1, :co_recalculated_a1, :??_acetone_def] => pv_calculation => :v1_acetone_area_corrected_fick)
    transform!(dataset, [:a2_acetone_area_volumecorrected, :g2_acetone_area_corrected, :ve_a2, :co_recalculated_a2, :??_acetone_def] => pv_calculation => :v2_acetone_area_corrected_fick)


    # Retentions

    retention(pa, pv) = @. pa / pv

    transform!(dataset, [:a1_sf6_area_volumecorrected, :v1_sf6_area_fick] => retention => :papv1_sf6)
    transform!(dataset, [:a2_sf6_area_volumecorrected, :v2_sf6_area_fick] => retention => :papv2_sf6)

    transform!(dataset, [:a1_ethane_area_volumecorrected, :v1_ethane_area_fick] => retention => :papv1_ethane)
    transform!(dataset, [:a2_ethane_area_volumecorrected, :v2_ethane_area_fick] => retention => :papv2_ethane)

    transform!(dataset, [:a1_cyclo_area_volumecorrected, :v1_cyclo_area_fick] => retention => :papv1_cyclo)
    transform!(dataset, [:a2_cyclo_area_volumecorrected, :v2_cyclo_area_fick] => retention => :papv2_cyclo)

    transform!(dataset, [:a1_iso_area_volumecorrected, :v1_iso_area_fick] => retention => :papv1_iso)
    transform!(dataset, [:a2_iso_area_volumecorrected, :v2_iso_area_fick] => retention => :papv2_iso)

    transform!(dataset, [:a1_dth_area_volumecorrected, :v1_dth_area_fick] => retention => :papv1_dth)
    transform!(dataset, [:a2_dth_area_volumecorrected, :v2_dth_area_fick] => retention => :papv2_dth)

    transform!(dataset, [:a1_acetone_area_volumecorrected, :v1_acetone_area_fick] => retention => :papv1_acetone)
    transform!(dataset, [:a2_acetone_area_volumecorrected, :v2_acetone_area_fick] => retention => :papv2_acetone)

    transform!(dataset, [:a1_acetone_area_volumecorrected, :v1_acetone_area_corrected_fick] => retention => :papv1_acetone_corrected)
    transform!(dataset, [:a2_acetone_area_volumecorrected, :v2_acetone_area_corrected_fick] => retention => :papv2_acetone_corrected)

    dataset.papv_sf6_ref .= 0.0001639
    dataset.papv_ethane_ref .= 0.034450
    dataset.papv_cyclo_ref .= 0.440427
    dataset.papv_iso_ref .= 0.750213
    dataset.papv_dth_ref .= 0.965825
    dataset.papv_acetone_ref .= 0.998811

    # Excretions

    retention(pe, pv) = @. pe / pv

    transform!(dataset, [:g1_sf6_area, :v1_sf6_area_fick] => retention => :pepv1_sf6)
    transform!(dataset, [:g2_sf6_area, :v2_sf6_area_fick] => retention => :pepv2_sf6)

    transform!(dataset, [:g1_ethane_area, :v1_ethane_area_fick] => retention => :pepv1_ethane)
    transform!(dataset, [:g2_ethane_area, :v2_ethane_area_fick] => retention => :pepv2_ethane)

    transform!(dataset, [:g1_cyclo_area, :v1_cyclo_area_fick] => retention => :pepv1_cyclo)
    transform!(dataset, [:g2_cyclo_area, :v2_cyclo_area_fick] => retention => :pepv2_cyclo)

    transform!(dataset, [:g1_iso_area, :v1_iso_area_fick] => retention => :pepv1_iso)
    transform!(dataset, [:g2_iso_area, :v2_iso_area_fick] => retention => :pepv2_iso)

    transform!(dataset, [:g1_dth_area, :v1_dth_area_fick] => retention => :pepv1_dth)
    transform!(dataset, [:g2_dth_area, :v2_dth_area_fick] => retention => :pepv2_dth)

    transform!(dataset, [:g1_acetone_area, :v1_acetone_area_fick] => retention => :pepv1_acetone)
    transform!(dataset, [:g2_acetone_area, :v2_acetone_area_fick] => retention => :pepv2_acetone)

    transform!(dataset, [:g1_acetone_area_corrected, :v1_acetone_area_corrected_fick] => retention => :pepv1_acetone_corrected)
    transform!(dataset, [:g2_acetone_area_corrected, :v2_acetone_area_corrected_fick] => retention => :pepv2_acetone_corrected)


    dataset.pepv_sf6_ref .= 0.0002
    dataset.pepv_ethane_ref .= 0.05
    dataset.pepv_cyclo_ref .= 0.20
    dataset.pepv_iso_ref .= 0.27
    dataset.pepv_dth_ref .= 0.34
    dataset.pepv_acetone_ref .= 0.37

    # Calculations for the partition coefficients

    dataset.vb_net = @. @with(dataset, (:vb - :c1_dry) / :pooled_blood_density)
    dataset.vg_net = dataset.vg .- dataset.vb_net

    dataset.vbx_net = @. @with(dataset, (:vbx - :syrx_dry) / :pooled_blood_density)
    dataset.vgx_net = dataset.vgx .- dataset.vbx_net

    partition_coefficient(vg, vb, s, sr) = @. vg / vb / (s / sr - 1)

    transform!(dataset, [:vg_net, :vb_net, :s_sf6_area, :sr_sf6_area] => partition_coefficient => :??_sf6)
    transform!(dataset, [:vg_net, :vb_net, :s_ethane_area, :sr_ethane_area] => partition_coefficient => :??_ethane)
    transform!(dataset, [:vg_net, :vb_net, :s_cyclo_area, :sr_cyclo_area] => partition_coefficient => :??_cyclo)
    transform!(dataset, [:vg_net, :vb_net, :s_iso_area, :sr_iso_area] => partition_coefficient => :??_iso)
    transform!(dataset, [:vg_net, :vb_net, :s_dth_area, :sr_dth_area] => partition_coefficient => :??_dth)
    transform!(dataset, [:vg_net, :vb_net, :s_acetone_area, :sr_acetone_area] => partition_coefficient => :??_acetone)


    return nothing


end

# Once the the output files of VQBOHR are reimported into the dataset, it is needed to choose the sample with the 
# smallest RSS. The next function takes the dataset as argument and returns a final dataset in which all the MIGET
# variables are also outputted with a _best suffix to indicate the ones with the lowest RSS. Morever also a list of
# pairs is outputted "tup" with a combination between names on the dataset and the associated best sample. This is useful 
# for later use. 

function select_best_sample(df :: DataFrame)   

    l = zeros(size(df, 1), 27)
    tup = []

    for (i,row) in enumerate(eachrow(df))
        
    if row.rss_a1 <= row.rss_a2
            row_df = DataFrame(row)
            df_best = row_df[:, Between(:bf_001_01_a1, :vent_dead_a1)]
            df_best.rss_best .= row.rss_a1
            m = values(df_best[1, :])
            tuple_to_add = (df.surname[i], "a1")
            push!(tup, tuple_to_add)
            
            
    elseif row.rss_a1 > row.rss_a2
            println("")
            row_df = DataFrame(row)
            df_best = row_df[:, Between(:bf_001_01_a2, :vent_dead_a2)]
            df_best.rss_best .= row.rss_a2
            m = values(df_best[1,:])
            tuple_to_add = (df.surname[i], "a2")
            push!(tup, tuple_to_add)

    end

    l[i, :] .= m

    end

    cols_keep = select(df, Between(:bf_001_01_a1, :vent_dead_a1))
    names_keep = names(cols_keep)
    new_names = replace.(names_keep, "_a1" => "_best")
    insert!(new_names, length(new_names)+1, "rss_best")

    df_merge = DataFrame(l, new_names, makeunique = true)
    select!(df_merge, Not(:rss_best_1))

    surs = df.surname
    df_merge[!, :surname] .= surs
    df_finale = outerjoin(df, df_merge, on = :surname)

    return df_finale, tup

end