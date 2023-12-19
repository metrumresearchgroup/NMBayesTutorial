#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2000/2000-3

/opt/NONMEM/nm75/run/nmfe75 2000-3.ctl  2000-3.lst  -parafile=2000-3.pnm -maxlim=2
