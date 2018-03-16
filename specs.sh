#!/bin/sh

# Create many PNG images and a Web page that links to them
# Open the file specs.html in a browser to view

path_to_files="/mnt/md0/total-power/2018"	# outrigger h5 files

echo -n > specs.html
for f in `ls $path_to_files/*.h5`
do 
  fout=`basename $f`
  echo "<p>$fout</p>" >> specs.html
  for pol in A B
  do
    if python 02_plot_spectra.py $f $pol
    then
      fout=`echo $f | sed s/\.h5/\_${pol}.png/`
      fout=`basename $fout`
      echo "<img src=spec_${fout}><br>" >> specs.html
    else echo "<p>$f $pol Failed</p>" >> specs.html
    fi
  done
  echo "<hr color=red>" >> specs.html
done
