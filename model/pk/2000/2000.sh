#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2000

/opt/NONMEM/nm75/run/nmfe75 2000.ctl  2000.lst  -maxlim=2
