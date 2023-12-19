#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2000pps

/opt/NONMEM/nm75/run/nmfe75 2000pps.ctl  2000pps.lst  -maxlim=2
