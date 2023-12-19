#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2000/2000-1

/opt/NONMEM/nm75/run/nmfe75 2000-1.ctl  2000-1.lst  -parafile=2000-1.pnm -maxlim=2
