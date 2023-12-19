#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2000/2000-2

/opt/NONMEM/nm75/run/nmfe75 2000-2.ctl  2000-2.lst  -parafile=2000-2.pnm -maxlim=2
