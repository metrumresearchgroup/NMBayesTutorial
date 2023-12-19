#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2005/2005-2

/opt/NONMEM/nm75/run/nmfe75 2005-2.ctl  2005-2.lst  -parafile=2005-2.pnm -maxlim=2
