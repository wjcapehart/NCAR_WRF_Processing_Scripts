#!/bin/csh
#
# LSF batch script to run the test MPI code
#
#BSUB -P USDM0004                       # Project number 99999999
#BSUB -a poe                            # select poe
#BSUB -n 1                            # number of total (MPI) tasks
#BSUB -R "span[ptile=1]"               # run a max of 32 tasks per node
#BSUB -J rtty_2050                            # job name
#BSUB -o rtty_2050_%J.out                      # output filename
#BSUB -e rtty_2050_%J.err                      # error filename
#BSUB -W 12:00                           # wallclock time
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
setenv YEAR 2050
module load ncl
ncl /glade/u/home/capehart/crunch_rtty_$YEAR.ncl >& /glade/u/home/capehart/rtty_$YEAR.txt

echo 
echo ===========================================
echo 
cat /glade/u/home/capehart/log_$YEAR.txt
 

echo
echo ===========================================
echo


echo
echo ===========================================
echo

echo "We are Out of Here Like Vladimir"
exit
