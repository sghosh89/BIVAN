#!/bin/bash
set -x
python subtract_SVL_and_to_nexus.py GA_Anolis_traits.csv  >trait_minus_svl.nexus 2>trait_minus_svl_contrasts.csv || exit
python pic1.py trait_minus_svl.nexus GA_Anolis_MCC.nex  >>trait_minus_svl_contrasts.csv || exit
python csv_to_nexus.py GA_Anolis_traits.csv > GA_Anolis_traits.nex 2>contrasts.csv || exit
python pic1.py GA_Anolis_traits.nex GA_Anolis_MCC.nex >> contrasts.csv || exit
