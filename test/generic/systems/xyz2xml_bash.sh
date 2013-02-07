#/bin/bash

# inputfile $1
# outputfile $2

if [ $# -lt 2 ]; then
	echo "Usage: $0 INPUTFILE OUTPUTFILE"
else
	echo '<?xml version="1.0" encoding="UTF-8"?>'               > $2
	echo "<!DOCTYPE scafacos_test SYSTEM 'scafacos_test.dtd'>" >> $2
	echo '<scafacos_test name="SYSTEM NAME" description="DESCRIPTION TEXT" reference_method="REFERENCE METHOD NAME" error_potential="1.0" error_field="1.0">' >> $2
	echo '  <configuration offset="0 0 0" box_a="0 0 0" box_b="0 0 0" box_c="0 0 0" periodicity="0 0 0" epsilon="metallic">' >> $2
	echo '    <duplicate times="1 1 1" rescale="0">' >> $2

	awk '{print "                <particle position=\""$1,$2,$3"\" q=\""$4"\"/>"}' < $1 >> $2

	echo '    </duplicate>'   >> $2
	echo '  </configuration>' >> $2
	echo '</scafacos_test>'   >> $2
fi

