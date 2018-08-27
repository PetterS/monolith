#  NLP written by GAMS Convert at 06/20/02 12:06:27
#  
#  Equation counts
#     Total       E       G       L       N       X
#        46      46       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        79      79       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#       327     113     214       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 50, >= 0, <= 1000;
var x3 := 50, >= 0, <= 1000;
var x4 := 50, >= 0, <= 1000;
var x5 := 50, >= 0, <= 1000;
var x6 := 50, >= 0, <= 1000;
var x7 := 50, >= 0, <= 1000;
var x8 := 50, >= 0, <= 1000;
var x9 := 50, >= 0, <= 1000;
var x10 := 50, >= 0, <= 1000;
var x11 := 50, >= 0, <= 1000;
var x12 := 0.2, >= 0, <= 10;
var x13 := 0.2, >= 0, <= 10;
var x14 := 0.2, >= 0, <= 10;
var x15 := 0.2, >= 0, <= 10;
var x16 := 0.2, >= 0, <= 10;
var x17 := 0.2, >= 0, <= 10;
var x18 := 0.2, >= 0, <= 10;
var x19 := 0.2, >= 0, <= 10;
var x20 := 0.2, >= 0, <= 10;
var x21 := 0.2, >= 0, <= 10;
var x22 := 100, >= 0, <= 1000;
var x23 := 100, >= 0, <= 1000;
var x24 := 100, >= 0, <= 1000;
var x25 := 100, >= 0, <= 1000;
var x26 := 100, >= 0, <= 1000;
var x27 := 0.2, >= 0, <= 10;
var x28 := 0.2, >= 0, <= 10;
var x29 := 0.2, >= 0, <= 10;
var x30 := 0.2, >= 0, <= 10;
var x31 := 0.2, >= 0, <= 10;
var x32 := 0.2, >= 0, <= 10;
var x33 := 0.2, >= 0, <= 10;
var x34 := 0.2, >= 0, <= 10;
var x35 := 0.2, >= 0, <= 10;
var x36 := 0.2, >= 0, <= 10;
var x37 := 50, >= 0, <= 1000;
var x38 := 50, >= 0, <= 1000;
var x39 := 50, >= 0, <= 1000;
var x40 := 50, >= 0, <= 1000;
var x41 := 50, >= 0, <= 1000;
var x42 := 50, >= 0, <= 1000;
var x43 := 50, >= 0, <= 1000;
var x44 := 50, >= 0, <= 1000;
var x45 := 50, >= 0, <= 1000;
var x46 := 50, >= 0, <= 1000;
var x47 := 50, >= 0, <= 1000;
var x48 := 50, >= 0, <= 1000;
var x49 := 50, >= 0, <= 1000;
var x50 := 50, >= 0, <= 1000;
var x51 := 50, >= 0, <= 1000;
var x52 := 50, >= 0, <= 1000;
var x53 := 50, >= 0, <= 1000;
var x54 := 50, >= 0, <= 1000;
var x55 := 50, >= 0, <= 1000;
var x56 := 50, >= 0, <= 1000;
var x57 := 50, >= 0, <= 1000;
var x58 := 50, >= 0, <= 1000;
var x59 := 50, >= 0, <= 1000;
var x60 := 50, >= 0, <= 1000;
var x61 := 50, >= 0, <= 1000;
var x62 := 50, >= 0, <= 1000;
var x63 := 50, >= 0, <= 1000;
var x64 := 50, >= 0, <= 1000;
var x65 := 50, >= 0, <= 1000;
var x66 := 50, >= 0, <= 1000;
var x67 := 100, >= 0, <= 1000;
var x68 := 0.2, >= 0, <= 10;
var x69 := 0.2, >= 0, <= 10;
var x70 := 1, >= 0, <= 10000;
var x71 := 1, >= 0, <= 10000;
var x72 := 1, >= 0, <= 10000;
var x73 := 1, >= 0, <= 10000;
var x74 := 1, >= 0, <= 10000;
var x75 >= 0, <= 10000;
var x76 >= 0, <= 10000;
var x77 >= 0, <= 10000;
var x78 >= 0, <= 10000;
var x79 >= 0, <= 10000;

minimize obj:  - x69;

subject to

e2:  - x2 - x3 - x4 - x5 - x6 = -100;

e3:  - x2 + x7 - x37 - x42 - x47 - x52 - x57 = 0;

e4:  - x3 + x8 - x38 - x43 - x48 - x53 - x58 = 0;

e5:  - x4 + x9 - x39 - x44 - x49 - x54 - x59 = 0;

e6:  - x5 + x10 - x40 - x45 - x50 - x55 - x60 = 0;

e7:  - x6 + x11 - x41 - x46 - x51 - x56 - x61 = 0;

e8: x12*x7 - (x27*x37 + x29*x42 + x31*x47 + x33*x52 + x35*x57) - 0.45*x2 = 0;

e9: x13*x7 - (x28*x37 + x30*x42 + x32*x47 + x34*x52 + x36*x57) - 0.55*x2 = 0;

e10: x14*x8 - (x27*x38 + x29*x43 + x31*x48 + x33*x53 + x35*x58) - 0.45*x3 = 0;

e11: x15*x8 - (x28*x38 + x30*x43 + x32*x48 + x34*x53 + x36*x58) - 0.55*x3 = 0;

e12: x16*x9 - (x27*x39 + x29*x44 + x31*x49 + x33*x54 + x35*x59) - 0.45*x4 = 0;

e13: x17*x9 - (x28*x39 + x30*x44 + x32*x49 + x34*x54 + x36*x59) - 0.55*x4 = 0;

e14: x18*x10 - (x27*x40 + x29*x45 + x31*x50 + x33*x55 + x35*x60) - 0.45*x5 = 0;

e15: x19*x10 - (x28*x40 + x30*x45 + x32*x50 + x34*x55 + x36*x60) - 0.55*x5 = 0;

e16: x20*x11 - (x27*x41 + x29*x46 + x31*x51 + x33*x56 + x35*x61) - 0.45*x6 = 0;

e17: x21*x11 - (x28*x41 + x30*x46 + x32*x51 + x34*x56 + x36*x61) - 0.55*x6 = 0;

e18:  - x7 + x22 = 0;

e19:  - x8 + x23 = 0;

e20:  - x9 + x24 = 0;

e21:  - x10 + x25 = 0;

e22:  - x11 + x26 = 0;

e23: x27*x22 - (x12*x7 - x70*x75) = 0;

e24: x28*x22 - (x13*x7 + x70*x75) = 0;

e25: x29*x23 - (x14*x8 - x71*x76) = 0;

e26: x30*x23 - (x15*x8 + x71*x76) = 0;

e27: x31*x24 - (x16*x9 - x72*x77) = 0;

e28: x32*x24 - (x17*x9 + x72*x77) = 0;

e29: x33*x25 - (x18*x10 - x73*x78) = 0;

e30: x34*x25 - (x19*x10 + x73*x78) = 0;

e31: x35*x26 - (x20*x11 - x74*x79) = 0;

e32: x36*x26 - (x21*x11 + x74*x79) = 0;

e33:  - x27*x28 + x75 = 0;

e34:  - x29*x30 + x76 = 0;

e35:  - x31*x32 + x77 = 0;

e36:  - x33*x34 + x78 = 0;

e37:  - x35*x36 + x79 = 0;

e38:    x22 - x37 - x38 - x39 - x40 - x41 - x62 = 0;

e39:    x23 - x42 - x43 - x44 - x45 - x46 - x63 = 0;

e40:    x24 - x47 - x48 - x49 - x50 - x51 - x64 = 0;

e41:    x25 - x52 - x53 - x54 - x55 - x56 - x65 = 0;

e42:    x26 - x57 - x58 - x59 - x60 - x61 - x66 = 0;

e43:  - x62 - x63 - x64 - x65 - x66 + x67 = 0;

e44: x67*x68 - (x62*x27 + x63*x29 + x64*x31 + x65*x33 + x66*x35) = 0;

e45: x67*x69 - (x62*x28 + x63*x30 + x64*x32 + x65*x34 + x66*x36) = 0;

e46:    x70 + x71 + x72 + x73 + x74 = 100;