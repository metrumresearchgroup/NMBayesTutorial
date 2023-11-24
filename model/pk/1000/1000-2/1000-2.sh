#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/1000/1000-2

/opt/NONMEM/nm75/run/nmfe75 1000-2.ctl  1000-2.lst  -parafile=1000-2.pnm -maxlim=2
