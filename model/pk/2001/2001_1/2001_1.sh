#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2001/2001_1

/opt/NONMEM/nm75/run/nmfe75 2001_1.ctl  2001_1.lst  -parafile=2001_1.pnm -maxlim=2
