#!/bin/bash
#########################################################################
# File Name    : generate-cnt.sh                                        #
# Author       : Neo                                                    #
# Mail         : niu.liu@nju.edu.cn                                     #
# Created Time : Wed 02 Jun 2021 06:58:31 PM CST                        #
#########################################################################

# Usage
help_msg(){
    echo "Usage: $0 <num_step> <cnt_file> <arc_file> <sou_file>"
    echo "Example: $0 20 opa2021s.cnt opa2021s.arc opa2021a.sou"
    echo "num_step : positive integer."
    echo "cnt_file : control file. I used the one on the OPAR machine."
    echo "           If you use a different one, you may need to modify line manually."
    echo "arc_file : arc file."
    echo "sou_file : .sou file outputed from Calc/Solve."
    echo "Note: this program should be run with generate-sou-list.py under the same directory."
}

# Check number of input argument
if [[ $# -ne 4 ]]
then
    help_msg
    exit 1
fi

nb_step=$1
cnt_file=$2
arc_file=$3
sou_file=$4

ipython generate-sou-list.py ${nb_step} ${sou_file}

for(( indx=1; indx<=nb_step; indx++ ))
do
    ind1=$((indx-1))
    sed -n '1,122p' ${cnt_file} > temp1
    cat list${ind1}  >> temp1
    sed -n '305,498p' ${cnt_file} >> temp1
    echo " ARCFILE                       opa2021s_${indx}.arc" >> temp1
    mv temp1 opa2021s_${indx}.cnt
    cp ${arc_file} opa2021s_${indx}.arc
    rm list${ind1}
done
