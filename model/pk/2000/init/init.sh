#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2000/init

/opt/NONMEM/nm75/run/nmfe75 init.ctl  init.lst  -parafile=init.pnm -maxlim=2
