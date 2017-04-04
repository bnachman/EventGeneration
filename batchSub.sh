#!/bin/bash

[ "$USER" == "pnef" ]     && WorkDir=/u/at/pnef/Work/Code/Reclustering/
[ "$USER" == "swiatlow" ] && WorkDir=/u/at/swiatlow/nfs/projects/Reclustering/
[ "$USER" == "bpn7" ]     && WorkDir=/nfs/slac/g/atlas/u01/users/bnachman/SLAC_pythia/Pileup
[ "$USER" == "jdamp" ]  && WorkDir= #Johannes fill in!
# add similar line if you are not pnef

SubFileLoc=`pwd`/_batchSingleSub.sh
#rm $SubFileLoc
DateSuffix=`date +%Y%m%d_%Hh%Mmin`

echo '#!/bin/bash
echo CD to $1
echo CMD is $2

cd $1
source /nfs/slac/g/atlas/u01/users/bnachman/SLAC_pythia/quick/setup.sh
cmd=$4

echo MAKING TEMP DIR $2
JOBFILEDIR=$2
mkdir $JOBFILEDIR
REALOUT=$3
echo MADE TEMP DIR $JOBFILEDIR
echo WILL COPY TO $REALOUT

shift
shift
echo Calling $cmd $*
$cmd $*
cp -r $JOBFILEDIR/*.txt $REALOUT
echo COPYING to $REALOUT
rm -rf $JOBFILEDIR
' > $SubFileLoc
chmod u+x $SubFileLoc

#----------------
Queue=long
nevents=100
njobs=100
npu_max=99
#180
#Queue=xlong
LogPrefix=`pwd`/logs/${DateSuffix}/${DateSuffix}_bsub_${mu}_
OutDirFinal=`pwd`/files/${DateSuffix}
mkdir -p `dirname $LogPrefix`
mkdir -p $OutDirFinal
echo
echo "Submitting $njobs jobs each with $nevents events to $Queue"
echo $LogPrefix
for (( npu=1; npu<=$npu_max; npu++ )) ;  do
    for (( ii=1; ii<=$njobs; ii++ )) ;  do
	echo $ii
	OutDir=/scratch/${DateSuffix}_${ii}/
	bsub -q ${Queue} -R 'select[(!preempt&&rhel60&&cvmfs&&inet)]' -o $LogPrefix${ii}.log $SubFileLoc           \
            ${WorkDir} ${OutDir} ${OutDirFinal} ./FCNC.exe  \
            --OutFile ${OutDir}/Sample_mu_${npu}_nevents_${nevents}_job_${ii}.txt \
            --NEvents ${nevents} --npu ${npu} 
    done
done