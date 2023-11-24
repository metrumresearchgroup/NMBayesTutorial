[ prob ]
1000

[ pkmodel ] cmt = "GUT,CENT,PERIPH", depot = TRUE

[ param ]
WT   = 70
ETA1 = 0
ETA2 = 0
ETA3 = 0
ETA4 = 0
ETA5 = 0

[ nmext ]
path = '../pk/1000/1000-1/1000-1.ext'
root = 'cppfile'

[ main ]
double V2WT   = log(WT/70);
double CLWT   = log(WT/70) * 0.75;
double V3WT   = log(WT/70);
double QWT    = log(WT/70) * 0.75;

capture KA = exp(THETA1 + ETA1 + ETA(1));
capture V2 = exp(THETA2 + V2WT + ETA2 + ETA(2));
capture CL = exp(THETA3 + CLWT + ETA3 + ETA(3));
capture V3 = exp(THETA4 + V3WT + ETA4 + ETA(4));
capture Q  = exp(THETA5 + QWT + ETA5 + ETA(5));

double S2 = V2/1000; //; dose in mcg, conc in mcg/mL

[ table ]
capture IPRED = CENT/S2;
capture Y = IPRED * (1 + EPS(1));
