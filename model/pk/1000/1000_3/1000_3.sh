#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1000/1000_3

/opt/NONMEM/nm75/run/nmfe75 1000_3.ctl  1000_3.lst  -parafile=1000_3.pnm -maxlim=2
