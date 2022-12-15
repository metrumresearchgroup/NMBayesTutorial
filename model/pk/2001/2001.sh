#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2001

/opt/NONMEM/nm75/run/nmfe75 2001.ctl  2001.lst  -maxlim=2
