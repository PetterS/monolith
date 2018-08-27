#  NLP written by GAMS Convert at 06/20/02 12:08:39
#  
#  Equation counts
#     Total       E       G       L       N       X
#        13      13       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        18      18       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        61      25      36       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 4.343494264, >= 4, <= 6;
var x2 := -4.313466584, >= -6, <= -4;
var x3 := 3.100750712, >= 2, <= 4;
var x4 := -2.397724192, >= -3, <= -1;
var x5 := 1.584424234, >= 1, <= 3;
var x6 := -1.551894266, >= -2, <= 0;
var x7 := 1.199661008, >= 0.5, <= 2.5;
var x8 := 0.212540694, >= -1.5, <= 0.5;
var x9 := 0.334227446, >= 0.2, <= 2.2;
var x10 := -0.199578662, >= -1.2, <= 0.8;
var x11 := 2.096235254, >= 0.1, <= 2.1;
var x12 := 0.057466756, >= -1.1, <= 0.9;
var x13 := 0.991133039, >= 0, <= 1;
var x14 := 0.762250467, >= 0, <= 1;
var x15 := 1.1261384966, >= 1.1, <= 1.3;
var x16 := 0.639718759, >= 0, <= 1;
var x17 := 0.159517864, >= 0, <= 1;

minimize obj: (x1 - 5)^2 + (5 + x2)^2 + (x3 - 3)^2 + (2 + x4)^2 + (x5 - 2)^2 + 
              (1 + x6)^2 + (x7 - 1.5)^2 + (0.5 + x8)^2 + (x9 - 1.2)^2 + (0.2 + 
              x10)^2 + (x11 - 1.1)^2 + (0.1 + x12)^2;

subject to

e2: x14/0.1570795**x15 - x1 + x13 = 0;

e3: x14/0.314159**x15 - x3 + x13 = 0;

e4: x14/0.4712385**x15 - x5 + x13 = 0;

e5: x14/0.628318**x15 - x7 + x13 = 0;

e6: x14/0.7853975**x15 - x9 + x13 = 0;

e7: x14/0.942477**x15 - x11 + x13 = 0;

e8:  - x17/0.1570795**x15 - x2 + 0.1570795*x16 = 0;

e9:  - x17/0.314159**x15 - x4 + 0.314159*x16 = 0;

e10:  - x17/0.4712385**x15 - x6 + 0.4712385*x16 = 0;

e11:  - x17/0.628318**x15 - x8 + 0.628318*x16 = 0;

e12:  - x17/0.7853975**x15 - x10 + 0.7853975*x16 = 0;

e13:  - x17/0.942477**x15 - x12 + 0.942477*x16 = 0;
