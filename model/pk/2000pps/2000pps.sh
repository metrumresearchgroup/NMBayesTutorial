#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2000pps

/opt/NONMEM/nm75/run/nmfe75 2000pps.ctl  2000pps.lst  -maxlim=2
