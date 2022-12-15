#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2000/2000_4

/opt/NONMEM/nm75/run/nmfe75 2000_4.ctl  2000_4.lst  -parafile=2000_4.pnm -maxlim=2
