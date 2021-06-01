#!/usr/local/bin/bash
#########################################################################
# File Name: replicate_arc.sh
# Author: Neo
# mail: niu.liu@nju.edu.cn
# Created Time: Tue May 25 15:20:15 2021
#########################################################################

OpaArc="../data/opa2021a.arc"
NjuArc="../data/nju_test.arc"
ComArc="../data/com_sess.arc"
NewArc="../data/new_sess.arc"
OutArc="../data/output.arc"

if [ -e ${OutArc}  ]; then
    rm ${OutArc}
fi

for sess in `cat ${ComArc}`
do
    grep ${sess} ${NjuArc} >> ${OutArc}
done

echo "** ----- SEPARATED LINE -----" >> ${OutArc}

for sess in `cat ${NewArc}`
do
    grep ${sess} ${OpaArc} >> ${OutArc}
done
