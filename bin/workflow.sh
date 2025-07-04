#!/bin/bash 

# a shell script to process a single EoS, solve for the macro file, extract features and (optionally) make plots in the SAME directory. 

# note: if outpath is specified anywhere, macropath etc will also have to be changed. 

EOS_PATH=$1           # path to the eos file    
MIN_PC2=$2            # min central pressure bound 
MAX_PC2=$3            # max central pressure bound

# filename for macro, defaults to eos-draw-{number}.csv-macro if not given via --outpath to integrate-tov

# 1. solve the tov equations 

./integrate-tov $EOS_PATH \
    $MIN_PC2 $MAX_PC2 \
    --central-eos-column baryon_density \
    --central-eos-column pressurec2 \
    --central-eos-column energy_densityc2 \
    --central-eos-column cs2c2 \
    --formalism logenthalpy \
    --interpolator-rtol 1e-8 \
    --integration-rtol 1e-8 \
    --V \
    --min-num-models 80

# if no default outpath given to integrate-tov, macro is written in the current directory as well 

MACROPATH="${EOS_PATH}-macro" 

# 2. extract features and (optionally) make a plot of the identified transitions

./extract-all-features $EOS_PATH \
    $MACROPATH \
    --verbose \
    --diff_arctan_threshold 0 \
    --cs_drop_threshold 0 \
    --diff_k_threshold 0.1 \
    --sw 2.5 \
    --plot

