#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2000/2000-4

/opt/NONMEM/nm75/run/nmfe75 2000-4.ctl  2000-4.lst  -parafile=2000-4.pnm -maxlim=2
