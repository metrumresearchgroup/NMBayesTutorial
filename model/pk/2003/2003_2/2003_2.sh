#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2003/2003_2

/opt/NONMEM/nm75/run/nmfe75 2003_2.ctl  2003_2.lst  -parafile=2003_2.pnm -maxlim=2
