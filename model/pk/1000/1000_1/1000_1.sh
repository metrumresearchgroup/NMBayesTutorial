#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1000/1000_1

/opt/NONMEM/nm75/run/nmfe75 1000_1.ctl  1000_1.lst  -parafile=1000_1.pnm -maxlim=2
