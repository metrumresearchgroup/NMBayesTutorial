#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2002/2002_1

/opt/NONMEM/nm75/run/nmfe75 2002_1.ctl  2002_1.lst  -parafile=2002_1.pnm -maxlim=2
