#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2000/2000_3

/opt/NONMEM/nm75/run/nmfe75 2000_3.ctl  2000_3.lst  -parafile=2000_3.pnm -maxlim=2
