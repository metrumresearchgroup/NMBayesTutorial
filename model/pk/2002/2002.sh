#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2002

/opt/NONMEM/nm75/run/nmfe75 2002.ctl  2002.lst  -maxlim=2
