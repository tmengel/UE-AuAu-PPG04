#!/bin/bash
currdir=/sphenix/user/tmengel/UE-AuAu-PPG04/offline/prelim



plots=( window/sigma_fit_cent0
        windowsys/sigma_fits_sys
        delta/area_mult_sub1_basic_0
        delta/area_mult_sub1_random_0
        rho/rho_comp
        delta/mult_basic_0
        delta/sub1_basic_0
        delta/area_basic_0
        RandomCones/sigma_et_vs_centrality
        RandomCones/sigma_et_vs_centrality_oversigma
        delta/area_mult_sub1_probe_0
        delta/area_mult_sub1_embed_0)

plotsName=( sPH-CONF-JET-2025-04-Fig1
            sPH-CONF-JET-2025-04-Fig2
            sPH-CONF-JET-2025-04-Fig3a
            sPH-CONF-JET-2025-04-Fig3b
            sPH-CONF-JET-2025-04-Fig4
            sPH-CONF-JET-2025-04-Fig5a
            sPH-CONF-JET-2025-04-Fig5b
            sPH-CONF-JET-2025-04-Fig5c
            sPH-CONF-JET-2025-04-Fig6 
            sPH-CONF-JET-2025-04-Fig7
            sPH-CONF-JET-2025-04-Fig8a
            sPH-CONF-JET-2025-04-Fig8b)

overLeafNames=(sigma_fit_cent0
               sigma_fits_sys
               area_mult_sub1_basic_0
               area_mult_sub1_random_0
               rho_comp
               mult_basic_0
               sub1_basic_0
               area_basic_0
               sigma_et_vs_centrality
               sigma_et_vs_centrality_oversigma
               area_mult_sub1_probe_0
               area_mult_sub1_embed_0)

auxplots=(  delta/area_basic_0_prelim
            delta/area_basic_0_prelim_randomized
            delta/area_basic_0_prelim_probe
            delta/area_basic_0_prelim_embed
            delta/mult_basic_0_prelim
            delta/mult_basic_0_prelim_randomized
            delta/mult_basic_0_prelim_probe
            delta/mult_basic_0_prelim_embed
            delta/sub1_basic_0_prelim
            delta/sub1_basic_0_prelim_randomized
            delta/sub1_basic_0_prelim_probe
            delta/sub1_basic_0_prelim_embed
            RandomCones/sigma_et_vs_centrality_Area
            RandomCones/sigma_et_vs_centrality_Multiplicity
            RandomCones/sigma_et_vs_centrality_Iterative)

auxNames=(  sPH-CONF-JET-2025-04-Aux1
            sPH-CONF-JET-2025-04-Aux2
            sPH-CONF-JET-2025-04-Aux3
            sPH-CONF-JET-2025-04-Aux4
            sPH-CONF-JET-2025-04-Aux5
            sPH-CONF-JET-2025-04-Aux6
            sPH-CONF-JET-2025-04-Aux7
            sPH-CONF-JET-2025-04-Aux8
            sPH-CONF-JET-2025-04-Aux9
            sPH-CONF-JET-2025-04-Aux10
            sPH-CONF-JET-2025-04-Aux11
            sPH-CONF-JET-2025-04-Aux12
            sPH-CONF-JET-2025-04-Aux13
            sPH-CONF-JET-2025-04-Aux14
            sPH-CONF-JET-2025-04-Aux15)



outputdir=$currdir/ppg04_final
mkdir -p "$outputdir"
mkdir -p "$outputdir/pdf" "$outputdir/png" "$outputdir/overleaf"
for i in "${!plots[@]}"; do
    
    plot=${currdir}/${plots[$i]}

    if [[ ! -f "$plot.pdf" ]]; then
        echo "Error: $plot.pdf does not exist. Skipping."
        continue
    fi
    outputfile="${outputdir}/pdf/${plotsName[$i]}.pdf"
    cp -f "$plot.pdf" "$outputfile"

    if [[ ! -f "$plot.png" ]]; then
        echo "Error: $plot.png does not exist. Skipping."
        continue
    fi
    outputfile_png="${outputdir}/png/${plotsName[$i]}.png"
    cp -f "$plot.png" "$outputfile_png"
    
    outputfile_ovrleaf="${outputdir}/overleaf/${overLeafNames[$i]}"
    cp -f "$plot.pdf" "${outputfile_ovrleaf}.pdf"
    cp -f "$plot.png" "${outputfile_ovrleaf}.png"

done

outputdir_aux=$currdir/ppg04_final/auxiliary
mkdir -p "$outputdir_aux"
mkdir -p "$outputdir_aux/pdf" "$outputdir_aux/png"

for i in "${!auxplots[@]}"; do
    
    auxplot=${currdir}/${auxplots[$i]}

    if [[ ! -f "$auxplot.pdf" ]]; then
        echo "Error: $auxplot.pdf does not exist. Skipping."
        continue
    fi
    outputfile_aux_pdf="${outputdir_aux}/pdf/${auxNames[$i]}.pdf"
    cp -f "$auxplot.pdf" "$outputfile_aux_pdf"

    if [[ ! -f "$auxplot.png" ]]; then
        echo "Error: $auxplot.png does not exist. Skipping."
        continue
    fi
    outputfile_aux_png="${outputdir_aux}/png/${auxNames[$i]}.png"
    cp -f "$auxplot.png" "$outputfile_aux_png"

done