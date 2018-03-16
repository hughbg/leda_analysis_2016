!#/bin/sh

# Run all the Python scripts. If you've already run script 00 and got the data then delete it.
for script in *.py
do 
  echo $script -------------------------
  if ! python $script data/outriggers_2016-01-27_14H10M25S.h5
  then exit
  fi
done
