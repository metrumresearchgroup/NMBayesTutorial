#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2003

/opt/NONMEM/nm75/run/nmfe75 2003.ctl  2003.lst  -maxlim=2
