#!/bin/sh


tmpfile=grid.txt

echo "" > $tmpfile

for x in 0.5 1.5 2.5 3.5; do
  for y in 0.5 1.5 2.5 3.5; do
    for z in 0.5 1.5 2.5 3.5; do
      echo "<particle position=\"$x $y $z\" q=\"1\"/>" >> $tmpfile
    done;
  done;
done;




