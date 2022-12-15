#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1000/1000_4

/opt/NONMEM/nm75/run/nmfe75 1000_4.ctl  1000_4.lst  -parafile=1000_4.pnm -maxlim=2
