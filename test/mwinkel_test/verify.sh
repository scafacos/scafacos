#!/bin/sh

RESULT=results

# setup-ids to be used (see setup.h)
SETUPS="0 1 2"
# Numbers of Particles to be used in tests
NPARTS="100 10000"
# (Pseudo)Periodicities to be used in tests
PERIODICITIES="000 111"

PROGNAMES="scafacos_verify scafacos_verify_f"

mkdir -p ${RESULT}

GREP="grep -s --binary-files=text"

for SETUP in ${SETUPS}
do
  for NPART in ${NPARTS}
  do
    for PERIODICITY in ${PERIODICITIES}
    do
      echo -e "\n############### NPARTICLES=${NPART}, PERIODICITY=${PERIODICITY}, SETUPID=${SETUP} ###############\n"
      echo -e "\t\t\t\tMETHOD\tSETUP\t\tEnergy\t\tTrace(Virial)\t\tSum(dot_product(r,E))"

      for i in `seq 0 7`
      do

        for PROGNAME in ${PROGNAMES}
        do
          FILENAME="${RESULT}/out.${PROGNAME}_${i}_${NPART}_${SETUP}_${PERIODICITY}.txt"
          rm -f ${FILENAME}
          COMMAND="${PROGNAME} ${i} ${NPART} ${SETUP} ${PERIODICITY}"
          ( ${COMMAND} ; true ) &> ${FILENAME}
          METHOD=`${GREP} "\-\->" ${FILENAME} | cut -d "-" -f 4 | tr -d ' '`
          SETUPA=`${GREP} "\=\=>" ${FILENAME} | cut -d "-" -f 2 | tr -d ' '`
          ENERGY=`${GREP} "energy/particle" ${FILENAME} | cut -d ":" -f 2 | tr -d ' '`
          TRAVIR=`${GREP} "trace of virial/particle" ${FILENAME} | cut -d ":" -f 2`
          TRAVI2=`${GREP} "dito (from forces)" ${FILENAME} | cut -d ":" -f 2`

          if [ ! -n "${ENERGY}" ]; then
            ENERGY="\t"`tac ${FILENAME} | head -n 3 | tr "\n" "\t" | cut -c 1-100`
          fi
           echo -e "${COMMAND}\t${METHOD}\t${SETUPA}\t\t${ENERGY}\t${TRAVIR}\t${TRAVI2}"
        done
      done
    done
  done
done

echo "############### Finished Solver Test ############################################################"
