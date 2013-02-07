#!/bin/sh -e

n=1000

make

./gendata grid_nacl_cube 512 > grid_nacl_cube_512.dat
./gendata grid_nacl_cube $n > grid_nacl_cube.dat
./gendata grid_nacl_cube 32768 > grid_nacl_cube_32K.dat

charge_type=0
./gendata halton_ellipsoid $n $charge_type > halton_ellipsoid_negpos.dat
./gendata halton_cylinder $n $charge_type > halton_cylinder_negpos.dat
./gendata hammersley_ball $n $charge_type > hammersley_ball_negpos.dat
./gendata hammersley_circle $n $charge_type > hammersley_circle_negpos.dat
./gendata hammersley_cube $n $charge_type > hammersley_cube_negpos.dat
./gendata hammersley_sphere $n $charge_type > hammersley_sphere_negpos.dat
./gendata hammersley_square $n $charge_type > hammersley_square_negpos.dat
./gendata plummer_ball $n $charge_type > plummer_ball_negpos.dat

# 1/64 = 0.015625
./gendata hammersley_two_balls $n $charge_type 0.015625 2 > hammersley_two_balls_negpos.dat
./gendata grid_face_centered_cube $n $charge_type > grid_face_centered_cube_negpos.dat
./gendata grid_body_centered_cube $n $charge_type > grid_body_centered_cube_negpos.dat

charge_type=1
./gendata halton_ellipsoid $n $charge_type > halton_ellipsoid_pos.dat
./gendata halton_cylinder $n $charge_type > halton_cylinder_pos.dat
./gendata hammersley_ball $n $charge_type > hammersley_ball_pos.dat
./gendata hammersley_circle $n $charge_type > hammersley_circle_pos.dat
./gendata hammersley_cube $n $charge_type > hammersley_cube_pos.dat
./gendata hammersley_sphere $n $charge_type > hammersley_sphere_pos.dat
./gendata hammersley_square $n $charge_type > hammersley_square_pos.dat
./gendata plummer_ball $n $charge_type > plummer_ball_pos.dat

# 1/64 = 0.015625
./gendata hammersley_two_balls $n $charge_type 0.015625 2 > hammersley_two_balls_pos.dat
./gendata grid_face_centered_cube $n $charge_type > grid_face_centered_cube_pos.dat
./gendata grid_body_centered_cube $n $charge_type > grid_body_centered_cube_pos.dat


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

	

