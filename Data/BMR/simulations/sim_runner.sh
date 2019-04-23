#!/bin/bash

pard="sims"
for stem in BivariateNormal ExtremeLeftTailDep ExtremeRightTailDep SomewhatLeftTailDep SomewhatRightTailDep 
do
    echo $stem
    if ! test -d "${pard}/${stem}" 
    then
        mkdir -p "${pard}/${stem}" || exit
    fi
    inp="MatsForMark/${stem}.csv"
    for ((rep=0; rep < 100; rep++))
    do
        outfn="${pard}/${stem}/rep${rep}.csv"
        errfn="${pard}/${stem}/err_${rep}_log.txt"
        echo "${stem} rep=${rep}"
        python phylo-trait-sim.py tree.nex "${inp}" >"${outfn}" 2>"${errfn}" || exit 
    done
done