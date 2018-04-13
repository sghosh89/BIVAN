#!/bin/bash
set -x
python subtract_SVL_and_to_nexus.py GA_Anolis_traits.csv  >trait_minus_svl.nexus 2>trait_minus_svl_contrasts.csv || exit
python pic1.py trait_minus_svl.nexus GA_Anolis_MCC.nex  >>trait_minus_svl_contrasts.csv || exit
