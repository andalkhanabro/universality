
#!/usr/bin/env bash

### A shell script to integrate each eos file FROM EOS_DIRECTORY and save the corresponding macrofiles in MACRO_DIRECTORY

source venv/bin/activate


EOS_DIRECTORY="eosfiles"          # directory relative to shell script where eosfiles are saved (i.e eosdata)
MACRO_DIRECTORY="macrofiles"      # directory relative to shell script where macro files would be saved (after integrating each to obtain the macro file via TOV solver)
EOSIDs_filename="eosIDlist.csv"   # csv file containing names of the csv "IDs" we want to work with.

MIN_PC2=1e13                      # min pressure we start integrating from
MAX_PC2=1e15                      # max pressure we integrate to

EOS_IDs=$(awk 'NR > 1 && NF {print $1}' $EOSIDs_filename | tr '\n' ' ') # to parse eos IDs chosen into a list

echo "IDs we'll generate macro files for:"  $EOS_IDs "\n" 

for id in $EOS_IDs; do
    EOS_FILENAME=$(printf "eos-draw-%06d.csv" "$id")
    EOS_FILEPATH="$EOS_DIRECTORY/$EOS_FILENAME"
    MACRO_FILENAME="$EOS_FILENAME-macro"
    MACRO_FILEPATH="$MACRO_DIRECTORY/$MACRO_FILENAME"

    echo ""
    echo $EOS_FILEPATH "\n"
    echo $MACRO_FILEPATH "\n"

    bin/integrate-tov $EOS_FILEPATH \
    $MIN_PC2 $MAX_PC2 \
    --central-eos-column baryon_density \
    --central-eos-column pressurec2 \
    --central-eos-column energy_densityc2 \
    --central-eos-column cs2c2 \
    --formalism logenthalpy \
    --interpolator-rtol 1e-8 \
    --integration-rtol 1e-8 \
    --V \
    --min-num-models 80 \
    --outpath $MACRO_FILEPATH    

done
