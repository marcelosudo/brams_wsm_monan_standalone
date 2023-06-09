#!/bin/csh

#---------------------------create the executable
rm -rf wsm.x
# run makefile
. mk.sh


#---------------------------create wsm.inp namelist
cat << Eof1 > wsm.inp

 &run
  prefix   ="check",
  mynum = 1,
  mcphys_type=5,
  time=3600,
 &end
Eof1

#--------------------------- run the job
./wsm.x >wsm.out




