#!/bin/sh -e

make

lengths="22 46 100"

for length in $lengths; do
  lengthcube=$(( length**3 ))
  file=grid_nacl_cube_periodic_$lengthcube.dat
  
  # generate geometry and charges
  ./gendata grid_nacl_cube_periodic $lengthcube > $file

  # convert .dat files into .xml format
  tail -n +3 $file > tmp_body.dat

  # convert to .xml and cut off header
  python ./nacl_xyz2xml.py tmp_body.dat | tail -n +6 > tmp_body.xml
  rm tmp_body.dat

  # create default header
  echo '<?xml version="1.0" ?>' > tmp_header.xml
  echo "<!DOCTYPE scafacos_test" >> tmp_header.xml
  echo "  SYSTEM 'scafacos_test.dtd'>" >> tmp_header.xml
  echo '<scafacos_test description="NaCl grid" error_field="1e-16" error_potential="1e-16" name="'$(basename $file .dat)'" reference_method="Analytic">' >> tmp_header.xml
  echo '	<configuration box_a="'$length'.0 0.0 0.0" box_b="0.0 '$length'.0 0.0" box_c="0.0 0.0 '$length'.0" epsilon="metallic" offset="0.0 0.0 0.0" periodicity="1 1 1">' >> tmp_header.xml

  # add header and zip the file
  cat tmp_header.xml tmp_body.xml > $(basename $file .dat).xml
  rm tmp_header.xml tmp_body.xml
  gzip $(basename $file .dat).xml
done

	

