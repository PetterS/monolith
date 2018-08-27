#  NLP written by GAMS Convert at 06/20/02 11:50:53
#  
#  Equation counts
#     Total       E       G       L       N       X
#        20      20       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        28      28       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        76      43      33       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 10, >= 10, <= 110;
var x2 := 10, >= 10, <= 110;
var x3 := 10, >= 10, <= 110;
var x4 := 10, >= 10, <= 110;
var x5 := 10, >= 10, <= 110;
var x6 := 10, >= 10, <= 110;
var x7 >= 0, <= 45;
var x8 >= 0, <= 45;
var x9 >= 0, <= 45;
var x10 >= 0, <= 45;
var x11 >= 0, <= 45;
var x12 >= 0, <= 45;
var x13 >= 0, <= 45;
var x14 >= 0, <= 45;
var x15 >= 0, <= 45;
var x16 >= 0, <= 45;
var x17 >= 0, <= 45;
var x18 >= 0, <= 45;
var x19 >= 0, <= 45;
var x20 >= 0, <= 45;
var x21 >= 0, <= 45;
var x22 := 100, >= 100, <= 200;
var x23 := 100, >= 100, <= 200;
var x24 := 100, >= 100, <= 200;
var x25 := 100, >= 100, <= 200;
var x26 := 100, >= 100, <= 200;
var x27 := 100, >= 100, <= 200;

minimize obj: 1300*(2000/(0.333333333333333*x1*x2 + 0.166666666666667*x1 + 
              0.166666666666667*x2))^0.6 + 1300*(1000/(0.666666666666667*x3*x4
               + 0.166666666666667*x3 + 0.166666666666667*x4))^0.6 + 1300*(1500
              /(0.666666666666667*x5*x6 + 0.166666666666667*x5 + 
              0.166666666666667*x6))^0.6;

subject to

e1:    x7 + x12 + x17 = 45;

e2:    x7 - x8 + x14 + x20 = 0;

e3:    x9 + x12 - x13 + x19 = 0;

e4:    x10 + x15 + x17 - x18 = 0;

e5:  - x8 + x9 + x10 + x11 = 0;

e6:  - x13 + x14 + x15 + x16 = 0;

e7:  - x18 + x19 + x20 + x21 = 0;

e8: x25*x14 + x27*x20 - x22*x8 + 100*x7 = 0;

e9: x23*x9 + x27*x19 - x24*x13 + 100*x12 = 0;

e10: x23*x10 + x25*x15 - x26*x18 + 100*x17 = 0;

e11: x8*x23 - x8*x22 = 2000;

e12: x13*x25 - x13*x24 = 1000;

e13: x18*x27 - x18*x26 = 1500;

e14:    x1 + x23 = 210;

e15:    x2 + x22 = 130;

e16:    x3 + x25 = 210;

e17:    x4 + x24 = 160;

e18:    x5 + x27 = 210;

e19:    x6 + x26 = 180;
