#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2003/2003_1

/opt/NONMEM/nm75/run/nmfe75 2003_1.ctl  2003_1.lst  -parafile=2003_1.pnm -maxlim=2
