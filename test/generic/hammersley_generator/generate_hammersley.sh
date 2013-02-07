#!/bin/sh -e

make

charge_type=1
./gendata hammersley_ball 10000 $charge_type > hammersley_ball_pos_1e4.dat
./gendata hammersley_ball 100000 $charge_type > hammersley_ball_pos_1e5.dat
./gendata hammersley_ball 1000000 $charge_type > hammersley_ball_pos_1e6.dat

# convert .dat files into .xml format
for file in *.dat; do
  tail -n +3 $file > tmp_body.dat

  # convert to .xml and cut off header
  python ../../generic/systems/xyz2xml.py tmp_body.dat | tail -n +6 > tmp_body.xml
  rm tmp_body.dat

  # create default header
  echo '<?xml version="1.0" ?>' > tmp_header.xml
  echo "<!DOCTYPE scafacos_test" >> tmp_header.xml
  echo "  SYSTEM 'scafacos_test.dtd'>" >> tmp_header.xml
  echo '<scafacos_test description="Short description of the testcase." error_field="1e-16" error_potential="1e-16" name="'$(basename $file .dat)'" reference_method="What method was used to generate the reference data.">' >> tmp_header.xml
  echo '	<configuration box_a="1.0 0.0 0.0" box_b="0.0 1.0 0.0" box_c="0.0 0.0 1.0" epsilon="metallic" offset="0.0 0.0 0.0" periodicity="0 0 0">' >> tmp_header.xml

  # add header and zip the file
  cat tmp_header.xml tmp_body.xml > $(basename $file .dat).noref.xml
  rm tmp_header.xml tmp_body.xml
  gzip $(basename $file .dat).noref.xml
done

	

