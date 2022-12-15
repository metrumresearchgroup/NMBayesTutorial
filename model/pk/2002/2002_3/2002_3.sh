#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2002/2002_3

/opt/NONMEM/nm75/run/nmfe75 2002_3.ctl  2002_3.lst  -parafile=2002_3.pnm -maxlim=2
