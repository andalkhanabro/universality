### A shell script to test batch processing logic for the new executable with tanh(K) functionality added to arctan(D) [old algorithm]

python3 bin/process2all-features eosIDlist.csv \
    --eos-dir outputs/my_test_data \
    --macro-dir outputs/my_test_data \
    --output-eos-dir outputs/my_test_output_andal \
    --eos-basename 'eos-draw-%(draw)06d.csv' \
    --macro-basename 'draw-macro-%(draw)06d.csv' \
    --macro2eos-central-baryon-density central_baryon_density baryon_density \
    --output-macro-column central_pressurec2 \
    --output-macro-column central_energy_densityc2 \
    --output-macro-column central_cs2c2 \
    -v

# each of these can be toggled to create directory structure for batch processing. however, to use the batch plotting functionality
# available at bin/plot_features.py, we need to modify the paths there that are used for dynamic construction and reading of plots. 


# CURRENTLY SET AT: 

#  --eos-dir outputs/my_test_data \
# --macro-dir outputs/my_test_data \
# --output-eos-dir outputs/my_test_output_andal \