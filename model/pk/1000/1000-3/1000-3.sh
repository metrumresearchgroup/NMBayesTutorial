#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/1000/1000-3

/opt/NONMEM/nm75/run/nmfe75 1000-3.ctl  1000-3.lst  -parafile=1000-3.pnm -maxlim=2
