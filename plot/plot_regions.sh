#!/bin/bash

cwd=`pwd`
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" #where this code is located
PARAMS_PLOT=params_plot.py
PLOT_STATIONS=$DIR/plot_stations.py #this code and plot_station.py are assumed to be in the same directory

opt_srf=""
opt_mparams=""

usage() { echo "Usage: $0 [-s srf_file] [-m model_params] [-c crns_file] params_plot_dir xyz_file" 1>&2; exit 1; }

while getopts ":hs:m:" arg; do
    case $arg in
        s) # -s OPTARG
            echo "-s $OPTARG"
            #check OPTARG is a valid file
            if [ -z "${OPTARG}" ] || [ ! -f `realpath ${OPTARG}` ]; then #if empty string or no such file
                usage
                exit 0
            else
                opt_srf="--srf $cwd/${OPTARG}"
            fi
            ;;
        m) # -m OPTARG
            #check OPTARG is a valid file
            echo "-m ${OPTARG}"
            if [ -z "${OPTARG}" ] || [ ! -f `realpath ${OPTARG}` ]; then #if empty string or no such file
                usage
                exit 0
            else
                opt_mparams="--model_params $cwd/${OPTARG}"
            fi
            ;;
		c)
			if [ -z "${OPTARG}" ] || [ ! -f `realpath ${OPTARG}` ]; then #if empty string or no such file
				usage
				exit 0
			else
				opt_cnrs="--srf_cnrs ${OPTARG}"
			fi
			;;
        h | *) #-h or everything else, show usage and exit
            usage 
            exit 0
            ;;
    esac
done

#parse positional arguments
regions_path=${@:$OPTIND:1} #first thing after processing all optional arguments
xyz=${@:$OPTIND+1:1} #the next thing

echo "================ Input argements ================"
echo "parmas_plot_dir: $regions_path"
echo "xyz_file: $xyz"
echo "srf option: $opt_srf"
echo "srf crns option: $opt_cnrs"
echo "model_params option: $opt_mparams"


regions_path=`realpath $regions_path` #get the real path even if relative path is given

for region_params_plot in $regions_path/params_plot*.py; do
    echo "================================================="
    echo "Getting ready to process $region_params_plot"
    rm $cwd/$PARAMS_PLOT*
    echo "Deleted existing $PARAMS_PLOT*"
    rm -rf GMT_WD_STATIONS
    ln -s $region_params_plot $cwd/$PARAMS_PLOT
    link_target=`readlink -f $cwd/$PARAMS_PLOT`
    echo "symlink $PARAMS_PLOT -> $link_target"

    echo "Executing: python $PLOT_STATIONS $xyz $opt_srf $opt_mparams $opt_cnrs"
    python $PLOT_STATIONS $xyz $opt_srf $opt_mparams $opt_cnrs

    output_img=`basename $region_params_plot` #remove PATH/ from PATH/params_plot_X.py
    output_img="${output_img#params_plot_}" #remove params_plot_ from params_plot_X.py
    output_img="${output_img%.*}" #remove .py from X.py
    i=0
    for old in $cwd/PNG_stations/c[0-9][0-9][0-9].png; do
        new=$output_img$(printf "%03d.png" "$i") #03 pad to length of 3
        echo "$old --> $new"
        new=$cwd/PNG_stations/$new
        mv $old $new
        let i=i+1
    done
done




