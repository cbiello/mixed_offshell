#!/bin/bash

runname="W_ATLAS"

## r0 = GamH / GammaHbb * BR(H->bb)
r0=1.0

# Clean old outputs
rm -f *avg*histo *summed*histo *dNNLO*histo

########################################
# Auto-detect scales
########################################
scalelist=$(ls ${runname}_nnlo_*_xMuR_*_xMuF_*_s*.histo 2>/dev/null \
    | sed -E 's/.*_xMuR_([0-9.]+)_xMuF_.*/\1/' | sort -u)

########################################
# Auto-detect channels
########################################
chlist=$(ls ${runname}_nnlo_*_xMuR_*_xMuF_*_s*.histo 2>/dev/null \
    | sed -E 's/.*_nnlo_([^_]+)_.*/\1/' | sort -u)

########################################
# Main loop
########################################
for scale in ${scalelist}; do
    echo ">>> Processing scale = ${scale}"

    for ch in ${chlist}; do
        echo "  -> Channel = ${ch}"

        ########################################
        # Auto-detect sections per (ch, scale)
        ########################################
        seclist=$(ls ${runname}_nnlo_${ch}_*_xMuR_${scale}_xMuF_${scale}_s*.histo 2>/dev/null \
            | sed -E 's/.*_nnlo_[^_]+_[^_]+_([^_]+)_[0-9]+_.*/\1/' | sort -u)

        for sec in ${seclist}; do
            echo "     * Section = ${sec}"

            myfiles=${runname}_nnlo_${ch}_*_${sec}_*_xMuR_${scale}_xMuF_${scale}_s*.histo

            avgfile=${runname}_avg_nnlo_${sec}_${ch}_xMuR_${scale}_xMuF_${scale}.histo
            r0file=${runname}_avg_nnlo_r0_${sec}_${ch}_xMuR_${scale}_xMuF_${scale}.histo

            ./rrsummer avg ${myfiles} ${avgfile}
            ./rrsummer mul $r0 ${avgfile} ${r0file}
        done

        ########################################
        # Sum over sections for this channel
        ########################################
        myfiles=${runname}_avg_nnlo_r0_*_${ch}_xMuR_${scale}_xMuF_${scale}.histo
        sumfiles=${runname}_summed_nnlo_${ch}_xMuR_${scale}_xMuF_${scale}.histo

        ./rrsummer add ${myfiles} ${sumfiles}
    done

    ########################################
    # Total NNLO correction (all channels)
    ########################################
    myfiles=${runname}_avg_nnlo_r0_*_*_xMuR_${scale}_xMuF_${scale}.histo
    sumfiles=${runname}_dNNLO_xMuR_${scale}_xMuF_${scale}.histo

    ./rrsummer add ${myfiles} ${sumfiles}

    ########################################
    # Add NLO → full NNLO
    ########################################
    if [ -f "${nlofile}" ]; then
        ./rrsummer add ${sumfiles} ${nlofile} ${nnlofile}
    else
        echo "WARNING: file not found, skipping NNLO combination"
    fi

done

########################################
# Optional: normalize bins
########################################
for scale in ${scalelist}; do
    sumfiles=${runname}_dNNLO_xMuR_${scale}_xMuF_${scale}.histo
    ./rrsummer sbn 1 ${sumfiles}
done

echo ">>> Done."




