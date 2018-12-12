#!/usr/bin/env bash



no_lines_to_filter=112

if [ -f thermo.csv ]; then

	rm thermo.csv

fi


tail -n+$no_lines_to_filter log.lammps>>thermo.csv
