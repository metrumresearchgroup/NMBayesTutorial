#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2004

/opt/NONMEM/nm75/run/nmfe75 2004.ctl  2004.lst  -maxlim=2
