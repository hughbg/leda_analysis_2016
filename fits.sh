#!/bin/sh

# Remove the fits/hkl files in the current directory and optionally
# in the h5 directory if you want to start over
doit() {
  path_to_files="$1"       # outrigger h5 files

  for f in `ls -t $path_to_files/*.h5`
  do 
    base=`echo $f | sed s/\.h5//`
    for ant in 252A 254A 254B 255A 255B
    do 
      if [ ! -r "${base}_${ant}.fits" ]
      then 
        echo ">>> Doing $f $ant fits"
        ./plot_ds9.py $2 $3 --median --no_show --flag $f $ant 
        mv *.fits $path_to_files
      else
        echo ">>> $f $ant done"
      fi
    done
    if [ ! -r "${base}.hkl" ]
    then
      echo ">>> Doing $f hickle"
      ./01_plot_waterfall.py $2 $3 --median --no_show --dump --flag $f 
      mv *.hkl $path_to_files
    else
      echo ">>> $f hkl done"
    fi
  done
}


rm /mnt/md0/total-power/2018/*.{fits,hkl}
doit /mnt/md0/total-power/2018 --new_cal
pushd /mnt/md0/total-power/2018
tar cf night.tar *.fits *.hkl
gzip night.tar
rm *.hkl *.fits
popd
rm *.fits *.hkl
doit /mnt/md0/total-power/2018 --new_cal --all_lsts
pushd /mnt/md0/total-power/2018
tar cf all_lsts.tar *.fits *.hkl
gzip all_lsts.tar
rm *.hkl *.fits
popd
rm *.fits *.hkl

rm /mnt/md0/total-power/*.{fits,hkl}
doit /mnt/md0/total-power
pushd /mnt/md0/total-power
tar cf night.tar *.fits *.hkl
gzip night.tar
rm *.hkl *.fits
popd
rm *.fits *.hkl
doit /mnt/md0/total-power --all_lsts
pushd /mnt/md0/total-power
tar cf all_lsts.tar *.fits *.hkl
gzip all_lsts.tar
rm *.hkl *.fits
popd
rm *.fits *.hkl

