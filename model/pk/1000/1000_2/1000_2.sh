#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1000/1000_2

/opt/NONMEM/nm75/run/nmfe75 1000_2.ctl  1000_2.lst  -parafile=1000_2.pnm -maxlim=2
