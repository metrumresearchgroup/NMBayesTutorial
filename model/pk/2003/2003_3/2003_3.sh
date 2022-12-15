#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2003/2003_3

/opt/NONMEM/nm75/run/nmfe75 2003_3.ctl  2003_3.lst  -parafile=2003_3.pnm -maxlim=2
