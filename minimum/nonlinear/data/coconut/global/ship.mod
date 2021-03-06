#  NLP written by GAMS Convert at 06/20/02 13:03:15
#  
#  Equation counts
#     Total       E       G       L       N       X
#        17       5      12       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        11      11       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        53      22      31       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 934.032;
var x2 := 934.032;
var x3 := 1011.868;
var x4 := 1.2, >= 1.05;
var x5 := 1.2, >= 1.05;
var x6 := 1.3, >= 1.05;
var x7 := 45.8;
var x8 := 43.2;
var x9 := 30.5;
var x10 := 76.3939536510076;

minimize obj: (0.0039*x7 + 0.0039*x8)*(495*x4 + 385*x5 + 315*x6)/x10;

subject to

e1:  - 0.5*x9*x4*(0.8*x7 + 0.333333333333333*x8) + x1 = 0;

e2:  - 0.5*x9*x5*(0.8*x7 + 0.333333333333333*x8) + x2 = 0;

e3:  - 0.5*x9*x6*(0.8*x7 + 0.333333333333333*x8) + x3 = 0;

e4: (x10 - x7)^2 - (x8^2 - x9^2) = 0;

e5:    x1 - 8.4652734375*x10 >= 0;

e6:    x2 - 9.65006510416667*x10 >= 0;

e7:    x3 - 8.8716796875*x10 >= 0;

e8: 0.5*x1*x9 - 2.2*(8.4652734375*x10)**1.33333333333333 >= 0;

e9: 0.5*x2*x9 - 2.2*(9.65006510416667*x10)**1.33333333333333 >= 0;

e10: 0.5*x3*x9 - 2.2*(8.8716796875*x10)**1.33333333333333 >= 0;

e11:    x4 - 0.0111771747883801*x7 >= 0.2;

e12:    x5 - 0.0137655360411427*x7 >= 0.2;

e13:    x6 - 0.0155663872253648*x7 >= 0.2;

e14:    x4 - 0.0111771747883801*x8 >= 0.2;

e15:    x5 - 0.0137655360411427*x8 >= 0.2;

e16:    x6 - 0.0155663872253648*x8 >= 0.2;
