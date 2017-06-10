#!/bin/csh
#
# LSF batch script to run the test MPI code
#
#BSUB -P USDM0004                       # Project number 99999999
#BSUB -a poe                            # select poe
#BSUB -n 1                            # number of total (MPI) tasks
#BSUB -R "span[ptile=1]"               # run a max of 32 tasks per node
#BSUB -J daily_rtty                            # job name
#BSUB -o daily_rtty_%J.out                      # output filename
#BSUB -e daily_rtty_%J.err                      # error filename
#BSUB -W 24:00                           # wallclock time
#BSUB -q geyser                        # queue
#BSUB -B
#BSUB -N



#
#source /usr/local/lsf/conf/cshrc.lsf
#
echo "Show Show Here we go!"
#
cd /glade/u/home/capehart/
#
#
source /glade/u/apps/opt/lmod/4.2.1/init/tcsh
#
setenv SSSS rtty 
module load ncl
ncl /glade/u/home/capehart/daily_${SSSS}.ncl >& /glade/u/home/capehart/daily_${SSSS}.txt

echo 
echo ===========================================
echo 
foreach YYYY (1990-2000 2020-2030 2030-2040 2050-2060 2080-2090)
nccopy -d 5 -k 3 -u -s ./daily_summary_${SSSS}_${YYYY}.nc  ./daily_summary_${SSSS}_${YYYY}.nc4
end

 

echo
echo ===========================================
echo


echo
echo ===========================================
echo

echo "We are Out of Here Like Vladimir"
exit
