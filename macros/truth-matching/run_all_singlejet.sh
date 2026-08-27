#!/bin/bash
version="08_23_2026_v001"
outdir="/sphenix/user/tmengel/dijet-ana-auau/macros/truth-matching/rootfiles/${version}/single_jet_flavor_fractions"
mkdir -p ${outdir}

for j in 10 20 30 ; do
    for type in scaled unscaled ; do 
        root -l -q -b "single_jet_flavor.C(\"/sphenix/user/tmengel/dijet-ana-auau/macros/truth-matching/rootfiles/${version}/jet${j}_hijing_${type}_all.root\", \"${outdir}/single_jet_flavor_fractions_jet${j}_${type}.root\", \"rho\")"
    done    
done
outdir="/sphenix/user/tmengel/dijet-ana-auau/macros/truth-matching/rootfiles/${version}/dijet_flavor_fractions"
mkdir -p ${outdir}
for j in 10 20 30 ; do
    for type in scaled unscaled ; do 
        root -l -q -b "leading_dijet_flavor.C(\"/sphenix/user/tmengel/dijet-ana-auau/macros/truth-matching/rootfiles/${version}/jet${j}_hijing_${type}_all.root\", \"${outdir}/leading_dijet_flavor_fractions_dijet${j}_${type}.root\", \"rho\")"
        done
done