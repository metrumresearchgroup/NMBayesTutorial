#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/1000

/opt/NONMEM/nm75/run/nmfe75 1000.ctl  1000.lst  -maxlim=2
