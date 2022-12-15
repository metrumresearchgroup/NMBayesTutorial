#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2000/2000_1

/opt/NONMEM/nm75/run/nmfe75 2000_1.ctl  2000_1.lst  -parafile=2000_1.pnm -maxlim=2
