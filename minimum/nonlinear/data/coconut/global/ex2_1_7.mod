#  NLP written by GAMS Convert at 06/20/02 11:41:50
#  
#  Equation counts
#     Total       E       G       L       N       X
#        11       1       0      10       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        21      21       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#       185     165      20       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0;
var x2 >= 0;
var x3 := 1.04289, >= 0;
var x4 >= 0;
var x5 >= 0;
var x6 >= 0;
var x7 >= 0;
var x8 >= 0;
var x9 >= 0;
var x10 >= 0;
var x11 := 1.74674, >= 0;
var x12 >= 0;
var x13 := 0.43147, >= 0;
var x14 >= 0;
var x15 >= 0;
var x16 := 4.43305, >= 0;
var x17 >= 0;
var x18 := 15.85893, >= 0;
var x19 >= 0;
var x20 := 16.4889, >= 0;

minimize obj:  - 0.5*((x1 - 2)^2 + 2*(x2 - 2)^2 + 3*(x3 - 2)^2 + 4*(x4 - 2)^2
               + 5*(x5 - 2)^2 + 6*(x6 - 2)^2 + 7*(x7 - 2)^2 + 8*(x8 - 2)^2 + 9*
              (x9 - 2)^2 + 10*(x10 - 2)^2 + 11*(x11 - 2)^2 + 12*(x12 - 2)^2 + 
              13*(x13 - 2)^2 + 14*(x14 - 2)^2 + 15*(x15 - 2)^2 + 16*(x16 - 2)^2
               + 17*(x17 - 2)^2 + 18*(x18 - 2)^2 + 19*(x19 - 2)^2 + 20*(x20 - 2
              )^2);

subject to

e2:  - 3*x1 + 7*x2 - 5*x4 + x5 + x6 + 2*x8 - x9 - x10 - 9*x11 + 3*x12 + 5*x13
     + x16 + 7*x17 - 7*x18 - 4*x19 - 6*x20 <= -5;

e3:    7*x1 - 5*x3 + x4 + x5 + 2*x7 - x8 - x9 - 9*x10 + 3*x11 + 5*x12 + x15
     + 7*x16 - 7*x17 - 4*x18 - 6*x19 - 3*x20 <= 2;

e4:  - 5*x2 + x3 + x4 + 2*x6 - x7 - x8 - 9*x9 + 3*x10 + 5*x11 + x14 + 7*x15
     - 7*x16 - 4*x17 - 6*x18 - 3*x19 + 7*x20 <= -1;

e5:  - 5*x1 + x2 + x3 + 2*x5 - x6 - x7 - 9*x8 + 3*x9 + 5*x10 + x13 + 7*x14
     - 7*x15 - 4*x16 - 6*x17 - 3*x18 + 7*x19 <= -3;

e6:    x1 + x2 + 2*x4 - x5 - x6 - 9*x7 + 3*x8 + 5*x9 + x12 + 7*x13 - 7*x14
     - 4*x15 - 6*x16 - 3*x17 + 7*x18 - 5*x20 <= 5;

e7:    x1 + 2*x3 - x4 - x5 - 9*x6 + 3*x7 + 5*x8 + x11 + 7*x12 - 7*x13 - 4*x14
     - 6*x15 - 3*x16 + 7*x17 - 5*x19 + x20 <= 4;

e8:    2*x2 - x3 - x4 - 9*x5 + 3*x6 + 5*x7 + x10 + 7*x11 - 7*x12 - 4*x13
     - 6*x14 - 3*x15 + 7*x16 - 5*x18 + x19 + x20 <= -1;

e9:    2*x1 - x2 - x3 - 9*x4 + 3*x5 + 5*x6 + x9 + 7*x10 - 7*x11 - 4*x12 - 6*x13
     - 3*x14 + 7*x15 - 5*x17 + x18 + x19 <= 0;

e10:  - x1 - x2 - 9*x3 + 3*x4 + 5*x5 + x8 + 7*x9 - 7*x10 - 4*x11 - 6*x12
      - 3*x13 + 7*x14 - 5*x16 + x17 + x18 + 2*x20 <= 9;

e11:    x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13
      + x14 + x15 + x16 + x17 + x18 + x19 + x20 <= 40;
