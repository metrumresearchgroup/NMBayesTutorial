$GLOBAL
#define cPRED (CENT/(V2)) * 1000 
  
[ prob ]
2000
  
[ param ]
FORM = 1
WT = 70

[ nmext ]
path = '../pk/2000/2000-1/2000-1.ext'
root = 'cppfile'

[ main ]
  
double LWT = log(WT/70); 
double CL  = exp(THETA1 + (LWT * 0.75) + ETA(1)); 
double V2  = exp(THETA2 + LWT + ETA(2));
double Q   = exp(THETA3 + (LWT * 0.75));
double V3  = exp(THETA4 + LWT);
double KA  = exp(THETA5);

[ pkmodel ]  cmt = "GUT CENT PERIPH", depot = TRUE
  
[ table ]
double IPRED = 1000 * (CENT / V2);
double Y = IPRED  * (1 + EPS(1));

[ capture ] cPRED IPRED Y KA V2 CL V3 Q

