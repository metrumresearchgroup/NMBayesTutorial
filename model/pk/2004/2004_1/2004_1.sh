#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2004/2004_1

/opt/NONMEM/nm75/run/nmfe75 2004_1.ctl  2004_1.lst  -parafile=2004_1.pnm -maxlim=2
