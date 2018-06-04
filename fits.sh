#!/bin/sh


path_to_files="/mnt/md0/total-power/2018"       # outrigger h5 files

for f in `ls -t $path_to_files/*.h5`
do 
  base=`echo $f | sed s/\.h5//`
  #for ant in 252A 252B 254A 254B 255A 255B
  #do 
  #  if [ ! -r "${base}_${ant}.fits ]
  #  then 
  #    echo ">>> Doing $f $ant fits"
  #    ./plot_ds9.py --median --no_show --flag $f $ant
  #    mv *.fits $path_to_files
  #  else
  #    echo ">>> $f $ant done"
  #  fi
  #done
  if [ ! -r "${base}.hkl" ]
  then
    echo ">>> Doing $f hickle"
    ./01_plot_waterfall.py --median --no_show --dump --flag $f
    mv *.hkl $path_to_files
  else
    echo ">>> $f hkl done"
  fi
done
