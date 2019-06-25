#! /bin/bash

if [ -z $1 ];then
	echo "./initiate.sh [initial-file-folder]"
	exit 0
fi

INIT=$1

mkdir -p topol tpr mdp gro analysis md_000

FN=`find $INIT/*.mdp`
tmp=`ls $FN | wc -l`
if [ $tmp -gt 1 ];then
	echo "[ERROR] I find more than 1 mdp files in your initial folder, exit."
	exit 0
fi
cp $FN ./mdp/md-nopbc.mdp

FN=`find $INIT/*.top`
tmp=`ls $FN | wc -l`
if [ $tmp -gt 1 ];then
	echo "[ERROR] I find more than 1 topol files in your initial folder, exit."
	exit 0
fi
cp $FN ./topol/topol_000.top

FN=`find $INIT/*.gro`
tmp=`ls $FN | wc -l`
if [ $tmp -gt 1 ];then
	echo "[ERROR] I found more than 1 gro files in your initial folder, exit."
	exit 0
fi
cp $FN ./md_000/out.gro

rm -rf `find .LASTIME_*`
touch .LASTIME_000
