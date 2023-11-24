$PROBLEM From bbr: see 1000-FOCE.yaml for details

$INPUT C NUM ID TIME SEQ CMT EVID AMT DV AGE WT HT EGFR ALB BMI SEX AAG
  SCR AST ALT CP TAFD TAD LDOS MDV BLQ PHASE

$DATA ../../../data/derived/analysis3.csv IGNORE=(C='C', BLQ=1)

$SUBROUTINE ADVAN4 TRANS4

$PK

;log transformed PK parms

  V2WT   = LOG(WT/70)
  CLWT   = LOG(WT/70) * 0.75
  V3WT   = LOG(WT/70)
  QWT    = LOG(WT/70) * 0.75

  MU_1   = THETA(1)
  MU_2   = THETA(2) + V2WT
  MU_3   = THETA(3) + CLWT
  MU_4   = THETA(4) + V3WT
  MU_5   = THETA(5) + QWT

  KA     = EXP(MU_1 + ETA(1))
  V2     = EXP(MU_2 + ETA(2))
  CL     = EXP(MU_3 + ETA(3))
  V3     = EXP(MU_4 + ETA(4))
  Q      = EXP(MU_5 + ETA(5))

  S2     = V2/1000 ; dose in mcg, conc in mcg/mL

$ERROR

  IPRED = F
  Y     = IPRED * (1 + EPS(1))

$THETA  ; log values
  (0.5)   ;  1 KA (1/hr) - 1.5
  (3.5)   ;  2 V2 (L) - 60
  (1)     ;  3 CL (L/hr) - 3.5
  (4)     ;  4 V3 (L) - 70
  (2)     ;  5 Q  (L/hr) - 4

$OMEGA BLOCK(3)
  0.2   ;ETA(KA)
  0.01 0.2   ;ETA(V2)
  0.01 0.01 0.2   ;ETA(CL)

$OMEGA
  0 FIX    ; ETA(V3)
  0 FIX    ; ETA(Q)

$SIGMA
  0.05     ; 1 pro error

$EST MAXEVAL=9999 METHOD=1 INTER SIGL=6 NSIG=3 PRINT=1 MSFO=./1000-FOCE.msf 
$COV PRINT=E
$TABLE NUM CL V2 Q V3 KA ETAS(1:LAST) IPRED NPDE CWRES NOPRINT ONEHEADER FILE=1000-FOCE.tab
