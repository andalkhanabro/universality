#!/bin/bash


# change these if your eos ID file is different, and your "flat" directories for eos and macro files are different

EOS_LIST="eosIDlist.csv"
ORIG_EOS_DIR="eosfiles"
ORIG_MAC_DIR="macrofiles" 

# where we will save this "new" test directory structure for batch processing 
TEST_DIR_BASE="outputs/my_test_data"

# using default from process2all-features, but should be whatever is set there. we dont have that many files rn. 

NUM_PER_DIR=1000

# refresh from any old runs 
echo "Cleaning up old test directory!"
rm -rf "$TEST_DIR_BASE"

# make new dir and build it 

mkdir -p "$TEST_DIR_BASE"
if [ ! -f "$EOS_LIST" ]; then
    echo "Error: EoS list '$EOS_LIST' not found!"
    exit 1
fi

echo "Building test directory structure in: $TEST_DIR_BASE"


tail -n +2 "$EOS_LIST" | while read -r eos_num; do
   

    eos_num=$(printf "%.0f" "$eos_num")

    moddraw=$((eos_num / NUM_PER_DIR))

    moddraw_str=$(printf "%06d" "$moddraw") # 6-digit padding

    subdir_name="DRAWmod${NUM_PER_DIR}-${moddraw_str}"
    target_dir="${TEST_DIR_BASE}/${subdir_name}"

    mkdir -p "$target_dir"

    eos_num_str=$(printf "%06d" "$eos_num")

    src_eos_file="${ORIG_EOS_DIR}/eos-draw-${eos_num_str}.csv"
    src_mac_file="${ORIG_MAC_DIR}/eos-draw-${eos_num_str}.csv-macro" 
    
    dest_eos_file="${target_dir}/eos-draw-${eos_num_str}.csv"
    dest_mac_file="${target_dir}/draw-macro-${eos_num_str}.csv"

    if [ -f "$src_eos_file" ] && [ -f "$src_mac_file" ]; then
        cp "$src_eos_file" "$dest_eos_file"
        cp "$src_mac_file" "$dest_mac_file" # copy and rename match csv convention? ask Reed abt this 
        echo "  Copied EoS $eos_num -> ${target_dir}/"
    else
        echo "  WARNING: Could not find source files for EoS $eos_num"
        echo "  (Checked for $src_eos_file and $src_mac_file)"
    fi
done

echo "built toy directory structure."