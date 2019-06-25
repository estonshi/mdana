#!/bin/bash

if [ -z $1 ];then
	if [ -z $LASTIME ];then
		echo "LASTIME is not defined, exit."
		exit 0
	fi
	TIMEP=$LASTIME
else
	TIMEP=$1
fi


GROMPP=$HOME/Documents/package/gromacs-4.6.7/bin/grompp_mpi
MDRUN=$HOME/Documents/package/gromacs-4.6.7/bin/mdrun_mpi


## generate new gro and top file
OLDGRO=`printf "md_%03d.gro" $TIMEP`
OLDTOP=`printf "topol_%03d.top" $TIMEP`
NEWTAG=`printf "%03d" $(($TIMEP+1))`
python edit_gro.py --gro ./gro/$OLDGRO --top ./topol/$OLDTOP --tag $NEWTAG
## generate new tpr file
$GROMPP -f ./mdp/md-nopbc.mdp -c ./gro/md_${NEWTAG}.ini.gro -p ./topol/topol_${NEWTAG}.top -o ./tpr/md_${NEWTAG}.tpr
if [ -f "./mdout.mdp" ];then
	rm mdout.mdp
fi
## copy files
mkdir md_$NEWTAG
cp ./tpr/md_${NEWTAG}.tpr ./md_$NEWTAG/md.tpr
cp ./submit.sh ./md_$NEWTAG/submit.sh
## run
cd md_$NEWTAG
qsub submit.sh
cd ..
