
LASTIMEF=`find .LASTIME_*`
if [ ! -n "$LASTIMEF" ];then
	echo "LASTIME is not defined, exit"
	exit 0
fi
LASTIME=`echo $LASTIMEF | tr -cd "[0-9]"`

JOB=`qstat -u ycshi`
OLDTAG=`printf "%03d" $LASTIME`

if [ -f "md_$OLDTAG/out.gro" ];then	
	cp ./md_$OLDTAG/out.gro ./gro/md_$OLDTAG.gro
	./load_job.sh $LASTIME
	LASTIME=$(($LASTIME+1))
	mv $LASTIMEF .LASTIME_$LASTIME
else
	echo "MD $LASTIME is not finished ..."
fi
		
