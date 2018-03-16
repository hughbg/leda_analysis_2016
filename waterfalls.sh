#!/bin/sh

# Create many PNG images and a Web page that links to them
# Open the file waterfalls.html in a browser to view

path_to_files="/mnt/md0/total-power/2018"       # outrigger h5 files

echo -n > waterfalls.html
for f in `ls $path_to_files/*.h5`
do 
  fout=`basename $f`
  echo "<p>$fout</p>" >> waterfalls.html
  for pol in A B
  do
    if python 01_plot_waterfall.py $f $pol
    then
      fout=`echo $f | sed s/\.h5/\_${pol}.png/`
      fout=`basename $fout`
      echo "<img src=$fout><br>" >> waterfalls.html
    else echo "<p>$f $pol Failed</p>" >> waterfalls.html
    fi
  done
  echo "<hr color=red>" >> waterfalls.html
done
