#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2002/2002_4

/opt/NONMEM/nm75/run/nmfe75 2002_4.ctl  2002_4.lst  -parafile=2002_4.pnm -maxlim=2
