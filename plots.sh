#!/bin/sh

# Create many PNG images and a Web page that links to them
# Open the file plots.html in a browser to view

path_to_files="/mnt/md0/total-power/2018"       # outrigger h5 files

echo -n > plots.html
for f in /mnt/md0/total-power/2018/outriggers_2018-03-10_11H56M37S.h5 #`ls $path_to_files/*.h5`
do 
  echo "$f --------------"
  fout=`basename $f`
  echo "<p>$fout</p>" >> plots.html
  for pol in A B
  do
    if python 01_plot_waterfall.py $f --pol $pol --save --noshow
    then
      fout=`echo $f | sed s/\.h5/\_${pol}.png/`
      fout=`basename $fout`
      echo "<img src=$fout><br>" >> plots.html
    else echo "<p>$f $pol waterfall Failed</p>" >> plots.html
    fi

    if python 01_plot_waterfall.py $f --pol $pol --flag --save --noshow
    then
      fout=`echo $f | sed s/\.h5/\_${pol}_flagged.png/`
      fout=`basename $fout`
      echo "Flagged<p><img src=$fout><br>" >> plots.html
    else echo "<p>$f $pol flagged waterfall Failed</p>" >> plots.html
    fi
  done

  if python plot_peaks.py $f --save --noshow
  then
    fout=`echo $f | sed s/\.h5/.png/`
    fout=`basename $fout`
    echo "<p>----------------<br>Max(channel) by channel<p><img src=peaks_$fout><br>" >> plots.html
  else echo "<p>$f $pol peaks Failed</p>" >> plots.html
  fi
  if python plot_peaks.py $f --save --noshow --flag
  then
    fout=`echo $f | sed s/\.h5/\_flagged.png/`
    fout=`basename $fout`
    echo "Max(channel) by channel Flagged<p><img src=peaks_$fout><br>" >> plots.html
  else echo "<p>$f $pol flagged peaks Failed</p>" >> plots.html
  fi
  echo "<hr color=red>" >> plots.html
done
exit


