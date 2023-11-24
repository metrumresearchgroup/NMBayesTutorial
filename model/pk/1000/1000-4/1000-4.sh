#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/1000/1000-4

/opt/NONMEM/nm75/run/nmfe75 1000-4.ctl  1000-4.lst  -parafile=1000-4.pnm -maxlim=2
