#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2000/2000_2

/opt/NONMEM/nm75/run/nmfe75 2000_2.ctl  2000_2.lst  -parafile=2000_2.pnm -maxlim=2
