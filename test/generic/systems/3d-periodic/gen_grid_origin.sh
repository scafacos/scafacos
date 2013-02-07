#!/bin/sh


tmpfile=grid_origin.txt

echo "" > $tmpfile

for x in 0 1 2 3; do
  for y in 0 1 2 3; do
    for z in 0 1 2 3; do
      echo "<particle position=\"$x $y $z\" q=\"1\"/>" >> $tmpfile
    done;
  done;
done;




