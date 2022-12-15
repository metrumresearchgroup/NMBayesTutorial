#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2002/2002_2

/opt/NONMEM/nm75/run/nmfe75 2002_2.ctl  2002_2.lst  -parafile=2002_2.pnm -maxlim=2
