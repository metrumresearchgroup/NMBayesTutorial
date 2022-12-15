#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2005

/opt/NONMEM/nm75/run/nmfe75 2005.ctl  2005.lst  -maxlim=2
