#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,amsfonts}'

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from universality.properties import features
from universality.utils import (io, utils)
import os
import traceback

########

INPUT_DATA_ROOT = "/Users/aandalkhanabro/Desktop/CITA 2025/universality/outputs/my_test_data"
OUTPUT_DATA_ROOT = "/Users/aandalkhanabro/Desktop/CITA 2025/universality/outputs/my_test_output_andal"
PLOT_SAVE_DIR = "/Users/aandalkhanabro/Desktop/CITA 2025/universality/outputs/plots"
EOS_NUM_PER_DIR = 1000 # NOTE: modify if changed in process2all-features

parser = ArgumentParser()


parser.add_argument(
    '--id-list', 
    type=str, 
    default='eosIDlist.csv', 
    help='csv file used in the batch processing'
)
parser.add_argument('-v', default=False, action='store_true')
parser.add_argument('-vd', default=False, action='store_true', help = "for deeper debugging (from feature extraction)")
parser.add_argument("--diff_k_threshold", default = features.DEFAULT_TANHK_DIFF_THRESHOLD, type=float, help = "min drop in K across the transition")
parser.add_argument("--cs_drop_threshold", default = features.DEFAULT_CS2C2_DROP_RATIO, type=float, help = "Rsc2 threshold")
parser.add_argument("--diff_arctan_threshold", default = features.DEFAULT_ARCTAN_DIFF_THRESHOLD, type=float, help = "min drop in arctan(D) across transition")
parser.add_argument("--sw", default = 3.5, type=float, help = "sigma for the gaussian smoothing applied to K")


args = parser.parse_args()
verbose = args.v


try:
    if verbose:
        print(f"reading EoS IDs from: {args.id_list}")
    #
    df = pd.read_csv(args.id_list)
    
    eos_id_list = df['EoS'].astype(int).tolist()
    
    if not eos_id_list:
        print(f"Error: No EoS IDs found in {args.id_list}. Exiting.")
        exit()

except FileNotFoundError:
    print(f"Error: ID list file not found at {args.id_list}. Exiting.")
    exit()
except KeyError:
    print(f"Error: CSV file {args.id_list} does not have an 'EoS' column. Exiting.")
    exit()
    

if verbose:
    print(f"Found {len(eos_id_list)} EoS IDs to plot.")
    print(f"plots will be saved to: {PLOT_SAVE_DIR}")



for eos_id in eos_id_list:
    
    if verbose:
        print(f"\n processing EoS id: {eos_id}")
    
    try:

        # dynamically construct paths, and predict what batch they will be in 

        batch_num = eos_id // EOS_NUM_PER_DIR
        batch_str = f"{batch_num:06d}"
        eos_id_str = f"{eos_id:06d}"

        # hardcoded (ask Reed to confirm logic) 

        batch_dir_name = f"DRAWmod{EOS_NUM_PER_DIR}-{batch_str}"
        eos_filename = f"eos-draw-{eos_id_str}.csv"
        macro_filename = f"draw-macro-{eos_id_str}.csv"
        features_filename = f"draw-macro-{eos_id_str}-all-features.csv"
        plot_filename = f"eos-draw-{eos_id_str}-plot.png"

        # eospaths and macropaths in the "data" directory, features in the "output" directory from batch processing for plotting 
        args.eospath = os.path.join(INPUT_DATA_ROOT, batch_dir_name, eos_filename)
        args.macropath = os.path.join(INPUT_DATA_ROOT, batch_dir_name, macro_filename)
        args.featurespath = os.path.join(OUTPUT_DATA_ROOT, batch_dir_name, features_filename)
        plotpath = os.path.join(PLOT_SAVE_DIR, plot_filename)
        
        if verbose:
            print(f" input path: {os.path.join(INPUT_DATA_ROOT, batch_dir_name)}")
            print(f" feature path: {args.featurespath}")

        
        # loading macro data
        macrodata, macro_columns = io.load(args.macropath, columns=['central_energy_densityc2', 'M', 'I', 'central_cs2c2', 
             'central_pressurec2', 'I', 'R', 'Lambda', 'central_baryon_density'])

        # loading eos data 
        eos_data, eos_cols = io.load(args.eospath, ["baryon_density", "cs2c2", "pressurec2"])
        baryon_density = eos_data[:,eos_cols.index("baryon_density")]
        eos_data = eos_data[np.argsort(baryon_density)] 
        eos_cs2c2 = eos_data[:,eos_cols.index("cs2c2")]
        eos_pressure = eos_data[:,eos_cols.index("pressurec2")]

        central_pressure = macrodata[:,macro_columns.index("central_pressurec2")]
        moi = macrodata[:,macro_columns.index("I")]
        mass = macrodata[:,macro_columns.index("M")]
        L = macrodata[:,macro_columns.index("Lambda")]
        macro_cs2c2 = macrodata[:,macro_columns.index("central_cs2c2")]
        central_energy_density = macrodata[:,macro_columns.index("central_energy_densityc2")]
        radius = macrodata[:,macro_columns.index("R")]
        rhoc = macrodata[:,macro_columns.index("central_baryon_density")]

        # to compute the arctan() transform for plotting 

        params, names, arctan_dlnR_dlnM, tanh_k = features.data2_tanh_features(
            rhoc, 
            mass=mass, 
            radius=radius, 
            eos_cs2c2=eos_cs2c2, 
            baryon_density=baryon_density,
            central_pressure=central_pressure,
            macro_data= macrodata, 
            eos_data=eos_data, 
            eos_cols = eos_cols, 
            macro_cols=macro_columns,
            diff_arctan_threshold=args.diff_arctan_threshold, 
            cs_drop_threshold=args.cs_drop_threshold, 
            diff_k_threshold=args.diff_k_threshold,
            tanh_k_sigma=args.sw,
            verbose=args.vd
        )

        # for plotting 


        #########

        mtov = np.argmax(mass)
        rhoc_at_mtov = rhoc[mtov]
        max_mass = mass[mtov] + 0.05

        fig, axes = plt.subplots(2, 2, figsize=(19, 15))
        ax1, ax2, ax3, ax4 = axes.flatten()  

        fullpath = args.featurespath
        extracted = pd.read_csv(fullpath)

        ### AX1 [eos space]
        ax1.plot(baryon_density, eos_cs2c2, color = "purple", lw=3)
        features.modify_ticks(ax1)
        ax1.set_xscale("log")

        start_pts_y = extracted["running_max_cs2c2_cs2c2"]
        start_pts_x = extracted["running_max_cs2c2_baryon_density"]

        end_pts_x = extracted["max_tanh_k_baryon_density"]
        end_pts_y = extracted["max_tanh_k_cs2c2"]

        ax1.plot(end_pts_x, end_pts_y, color= "blue", markerfacecolor='none', 
                                markeredgecolor="blue", markersize=13, linestyle='None', markeredgewidth=4, marker = 4)
        ax1.plot(start_pts_x, start_pts_y, color= "blue", markerfacecolor='none', 
                                markeredgecolor="red", markersize=13, linestyle='None', markeredgewidth=4, marker = 5)

        # annotate the identified transitions 

        for i in range(len(start_pts_x)):
            
            current_start_x = start_pts_x.iloc[i]
            current_end_x = end_pts_x.iloc[i]
            
            ax1.axvspan(min(current_start_x, current_end_x), 
                        max(current_start_x, current_end_x), 
                        color='orange', alpha=0.15, zorder=-1) 
            ax1.axvline(current_start_x, color='darkorange', linewidth=1.2, zorder=0)
            ax1.axvline(current_end_x, color='darkorange', linewidth=1.2, zorder=0)

        ax1.scatter(rhoc[mtov], macro_cs2c2[mtov], color = "black", alpha = 0.3, s= 200)
        macro_cs2c2_at_mtov = macro_cs2c2[mtov]
        rhoc_at_bound = rhoc[mtov]
        macro_cs2c2_at_bound = macro_cs2c2[mtov]
        epsilon = 1e15
        ax1.set_xlim(8e13, rhoc_at_bound + epsilon)
        ax1.scatter(rhoc_at_mtov, macro_cs2c2_at_mtov, color = "black", alpha = 0.5, s= 200)
        ax1.set_xlabel(r"$\rho$")
        ax1.set_ylabel(r"$c_s^2/c^2$")
        

        ### AX1 ENDS HERE 

        ## AX2 [mass vs arctan (old feature space)]
        ax2.plot(mass, arctan_dlnR_dlnM, color = "purple", lw = 3)
        features.modify_ticks(ax2)

        start_pts_x = extracted["running_max_cs2c2_M"]
        start_pts_y = extracted["running_max_cs2c2_arctan_dlnR_dlnM"]

        end_pts_x = extracted["max_tanh_k_M"]
        end_pts_y = extracted["max_tanh_k_arctan_dlnR_dlnM"]

        for i in range(len(start_pts_x)):
            
            current_start_x = start_pts_x.iloc[i]
            current_end_x = end_pts_x.iloc[i]
            
            ax2.axvspan(min(current_start_x, current_end_x), 
                        max(current_start_x, current_end_x), 
                        color='orange', alpha=0.15, zorder=-1) 
            ax2.axvline(current_start_x, color='darkorange', linewidth=1.2, zorder=0)
            ax2.axvline(current_end_x, color='darkorange', linewidth=1.2, zorder=0)

        ax2.plot(end_pts_x, end_pts_y, color= "blue", markerfacecolor='none', 
                                markeredgecolor="blue", markersize=13, linestyle='None', markeredgewidth=4, marker = 4)
        ax2.plot(start_pts_x, start_pts_y, color= "blue", markerfacecolor='none', 
                                markeredgecolor="red", markersize=13, linestyle='None', markeredgewidth=4, marker = 5)
        ax2.scatter(mass[mtov], arctan_dlnR_dlnM[mtov], color = "black", alpha = 0.5, s= 200)

        y_min, y_max = ax2.get_ylim()
        ax2.set_ylim(y_min, y_max) 
        ax2.axhspan(np.pi/2, y_max, alpha=0.5, hatch = "\\", color = "gray")
        ax2.axhspan(y_min, -np.pi/2, alpha=0.5, hatch = "\\", color = "gray")
        ax2.set_xlabel(r"$M_{\odot}$")
        ax2.set_ylabel(r"$\arctan(D)$")

        # AX2 ENDS HERE 

        ## AX3 [rho vs sound-speed squared]

        start_pts_x = extracted["running_max_cs2c2_baryon_density"]
        start_pts_y = extracted["running_max_cs2c2_pressurec2"]

        end_pts_x = extracted["max_tanh_k_baryon_density"]
        end_pts_y = extracted["max_tanh_k_pressurec2"]

        for i in range(len(start_pts_x)):
            
            current_start_x = start_pts_x.iloc[i]
            current_end_x = end_pts_x.iloc[i]
            
            ax3.axvspan(min(current_start_x, current_end_x), 
                        max(current_start_x, current_end_x), 
                        color='orange', alpha=0.15, zorder=-1) 
            ax3.axvline(current_start_x, color='darkorange', linewidth=1.2, zorder=0)
            ax3.axvline(current_end_x, color='darkorange', linewidth=1.2, zorder=0)

        ax3.plot(end_pts_x, end_pts_y, color= "blue", markerfacecolor='none', 
                                markeredgecolor="blue", markersize=13, linestyle='None', markeredgewidth=4, marker = 4)
        ax3.plot(start_pts_x, start_pts_y, color= "blue", markerfacecolor='none', 
                                markeredgecolor="red", markersize=13, linestyle='None', markeredgewidth=4, marker = 5)

        rhoc_at_mtov = rhoc[mtov]
        pressure_at_mtov = central_pressure[mtov] 

        ax3.plot(baryon_density[baryon_density <= rhoc_at_mtov], eos_pressure[eos_pressure <= pressure_at_mtov], color = "purple", lw = 3)
        ax3.scatter(baryon_density[baryon_density <= rhoc_at_mtov][-1], eos_pressure[eos_pressure <= pressure_at_mtov][-1], color = "black", alpha = 0.3, \
                    s = 200)

        epsilon = 1e14
        features.modify_ticks(ax3)
        ax3.set_xlim(1e14, baryon_density[baryon_density <= rhoc_at_mtov][-1] + epsilon)
        features.modify_ticks(ax3)
        ax3.set_xscale("log")
        ax3.set_yscale("log")
        ax3.set_ylim(1e12, 1e16)
        ax3.set_xlabel(r"$\rho$")
        ax3.set_ylabel(r"$p/c^2$")

        ## AX3 ENDS HERE 

        ### AX4 [k-space [RADIUS] vs central pressure]
        start_x, start_y = extracted["running_max_cs2c2_central_pressurec2"], extracted["running_max_cs2c2_tanh_k"]
        end_x, end_y = extracted["max_tanh_k_central_pressurec2"], extracted["max_tanh_k_tanh_k"]

        for i in range(len(start_x)):
            
            current_start_x = start_x.iloc[i]
            current_end_x = end_x.iloc[i]
            
            ax4.axvspan(min(current_start_x, current_end_x), 
                        max(current_start_x, current_end_x), 
                        color='orange', alpha=0.15, zorder=-1) 
            ax4.axvline(current_start_x, color='darkorange', linewidth=1.2, zorder=0)
            ax4.axvline(current_end_x, color='darkorange', linewidth=1.2, zorder=0)

        epsilon_for_plotting = 1e14
        ax4.plot(central_pressure[0:mtov] , tanh_k[0:mtov], lw = 3, color = "purple")
        features.modify_ticks(ax4)
        ax4.set_xscale("log")
        ax4.scatter(central_pressure[mtov], tanh_k[mtov], color = "black", alpha = 0.5, s= 200)
        ax4.plot(end_x, end_y, color= "blue", markerfacecolor='none', 
                                markeredgecolor="blue", markersize=13, linestyle='None', markeredgewidth=4, marker = 4)
        ax4.plot(start_x, start_y, color= "blue", markerfacecolor='none', 
                                markeredgecolor="red", markersize=13, linestyle='None', markeredgewidth=4, marker = 5)

        y_min, y_max = ax4.get_ylim()
        ax4.set_ylim(y_min, y_max) 
        ax4.axhspan(0, y_max, alpha=0.25, hatch = ".", color = "lightblue")
        ax4.axhspan(y_min, 0, alpha=0.25, hatch = ".", color = "lightcoral")
        ax4.set_xlabel(r"$P_c$")
        ax4.set_ylabel(r"$tanh(K_R)$")

        ### AX4 ENDS HERE 

        plt.tight_layout()
        plt.savefig(plotpath)

        if verbose:
            print(f" plot for {eos_id} saved to {plotpath}")

    except FileNotFoundError:
        print(f"  [!] WARNING: could not find all files for the specific EoS {eos_id}. skipping.")
        if verbose:
            # Print full error stack trace for debugging file path issues
            traceback.print_exc()
            
    except Exception as e:
        print(f"   error: failed to plot EoS {eos_id} due to an unexpected error?: {e}")
        if verbose:
        
            traceback.print_exc()

    finally:
        plt.close('all')

print(f"\n All {len(eos_id_list)} algorithm results have been plotted.")