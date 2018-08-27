rm hack.bas
/bin/umodem hack.bas STOP
V        REM  +------------------------------------------------+
X        REM  | HACK.BAS      (c) 19100   fr33 v4r14bl3z       |
XV       REM  |                                                |
XX       REM  | Brute-forces passwords on UM vIX.0 systems.    |
XXV      REM  | Compile with Qvickbasic VII.0 or later:        |
XXX      REM  |    /bin/qbasic hack.bas                        |
XXXV     REM  | Then run:                                      |
XL       REM  |   ./hack.exe username                          |
XLV      REM  |                                                |
L        REM  | This program is for educational purposes only! |
LV       REM  +------------------------------------------------+
LX       REM
LXV      IF ARGS() > I THEN GOTO LXXXV
LXX      PRINT "usage: ./hack.exe username"
LXXV     PRINT CHR(X)
LXXX     END
LXXXV    REM
XC       REM  get username from command line
XCV      DIM username AS STRING
C        username = ARG(II)
CV       REM  common words used in passwords
CX       DIM pwdcount AS INTEGER
CXV      pwdcount = LIII
CXX      DIM words(pwdcount) AS STRING
CXXV     words(I) = "airplane"
CXXX     words(II) = "alphabet"
CXXXV    words(III) = "aviator"
CXL      words(IV) = "bidirectional"
CXLV     words(V) = "changeme"
CL       words(VI) = "creosote"
CLV      words(VII) = "cyclone"
CLX      words(VIII) = "december"
CLXV     words(IX) = "dolphin"
CLXX     words(X) = "elephant"
CLXXV    words(XI) = "ersatz"
CLXXX    words(XII) = "falderal"
CLXXXV   words(XIII) = "functional"
CXC      words(XIV) = "future"
CXCV     words(XV) = "guitar"
CC       words(XVI) = "gymnast"
CCV      words(XVII) = "hello"
CCX      words(XVIII) = "imbroglio"
CCXV     words(XIX) = "january"
CCXX     words(XX) = "joshua"
CCXXV    words(XXI) = "kernel"
CCXXX    words(XXII) = "kingfish"
CCXXXV   words(XXIII) = "(\b.bb)(\v.vv)"
CCXL     words(XXIV) = "millennium"
CCXLV    words(XXV) = "monday"
CCL      words(XXVI) = "nemesis"
CCLV     words(XXVII) = "oatmeal"
CCLX     words(XXVIII) = "october"
CCLXV    words(XXIX) = "paladin"
CCLXX    words(XXX) = "pass"
CCLXXV   words(XXXI) = "password"
CCLXXX   words(XXXII) = "penguin"
CCLXXXV  words(XXXIII) = "polynomial"
CCXC     words(XXXIV) = "popcorn"
CCXCV    words(XXXV) = "qwerty"
CCC      words(XXXVI) = "sailor"
CCCV     words(XXXVII) = "swordfish"
CCCX     words(XXXVIII) = "symmetry"
CCCXV    words(XXXIX) = "system"
CCCXX    words(XL) = "tattoo"
CCCXXV   words(XLI) = "thursday"
CCCXXX   words(XLII) = "tinman"
CCCXXXV  words(XLIII) = "topography"
CCCXL    words(XLIV) = "unicorn"
CCCXLV   words(XLV) = "vader"
CCCL     words(XLVI) = "vampire"
CCCLV    words(XLVII) = "viper"
CCCLX    words(XLVIII) = "warez"
CCCLXV   words(XLIX) = "xanadu"
CCCLXX   words(L) = "xyzzy"
CCCLXXV  words(LI) = "zephyr"
CCCLXXX  words(LII) = "zeppelin"
CCCLXXXV words(LIII) = "zxcvbnm"
CCCXC    REM try each password
CCCXCV   PRINT "attempting hack with " + pwdcount + " passwords " + CHR(X)
CD       DIM i AS INTEGER
CDV      i = I
CDX      IF CHECKPASS(username, words(i)) THEN GOTO CDXXX
CDXV     i = i + I
CDXX     IF i > pwdcount THEN GOTO CDXLV
CDXXV    GOTO CDX
CDXXX    PRINT "found match!! for user " + username + CHR(X)
CDXXXV   PRINT "password: " + words(i) + CHR(X)
CDXL     END
CDXLV    PRINT "no simple matches for user " + username + CHR(X)
CDL      REM
CDLV     REM  the above code will probably crack passwords for many
CDLX     REM  users so I always try it first. when it fails, I try the
CDLXV    REM  more expensive method below.
CDLXX    REM
CDLXXV   REM  passwords often take the form
CDLXXX   REM    dictwordDD
CDLXXXV  REM  where DD is a two-digit decimal number. try these next:
MI        i = I
MII       DIM j AS INTEGER
MV        j = I
MVI       REM DIM k as INTEGER
MVII      REM k = I
MX        DIM pass AS STRING
MXV       pass = words(i) + CHR(XLVII + j) + CHR(XLVII + j) 
MXVI      PRINT "Trying " + pass + CHR(X)
MXX       IF CHECKPASS(username, pass) THEN GOTO MLXV 
MXXI      REM k = k + I  
MXXII     REM IF k > X THEN GOTO MXXV  
MXXIII    REM GOTO MX
MXXV      j = j + I  
MXXX      IF j > X THEN GOTO MXL  
MXXXV     GOTO MVI  
MXL       i = i + I  
MXLV      IF i > pwdcount THEN GOTO MXC  
ML        GOTO MII
MLV       REM  
MLX       REM
MLXV      PRINT "found advanced match!! for user " + username + CHR(X) 
MLXX      PRINT "password: " + pass + CHR(X) 
MLXXV     REM 
MLXXX     REM 
MLXXXV    REM 
MXC       PRINT "no advanced matches for user " + username + CHR(X) 
MXCV      REM 
MC        REM 
MCV       REM 
MCX       REM 
MCXV      REM 
MCXX      REM 
MCXXV     REM 
MCXXX     REM 
MCXXXV    REM 
MCXL      REM 
MCXLV     REM 
MCL       REM 
MCLV      REM 
MCLX      REM 
MCLXV     REM 
MCLXX     REM 
MCLXXV    REM 
MCLXXX    REM 
MCLXXXV   REM 
MCXC      REM 
MCXCV     REM 
MCC       REM 
MCCV      REM 
MCCX      REM 
MCCXV     REM 
MCCXX     REM 
MCCXXV    REM 
MCCXXX    REM 
MCCXXXV   REM 
MCCXL     REM 
MCCXLV    REM 
MCCL      REM 
MCCLV     REM 
MCCLX     REM 
MCCLXV    REM 
MCCLXX    REM 
MCCLXXV   REM 
MCCLXXX   REM 
MCCLXXXV  REM 
MCCXC     REM 
MCCXCV    REM 
MCCC      REM 
MCCCV     REM 
MCCCX     REM 
MCCCXV    REM 
MCCCXX    REM 
MCCCXXV   REM 
MCCCXXX   REM 
MCCCXXXV  REM 
MCCCXL    REM 
MCCCXLV   REM 
MCCCL     REM 
MCCCLV    REM 
MCCCLX    REM 
MCCCLXV   REM 
MCCCLXX   REM 
MCCCLXXV  REM 
MCCCLXXX  REM 
MCCCLXXXV REM 
MCCCXC    REM 
MCCCXCV   REM 
MCD       REM 
MCDV      REM 
MCDX      REM 
MCDXV     REM 
MCDXX     REM 
MCDXXV    REM 
MCDXXX    REM 
MCDXXXV   REM 
MCDXL     REM 
MCDXLV    REM 
MCDL      REM 
MCDLV     REM 
MCDLX     REM 
MCDLXV    REM 
MCDLXX    REM 
MCDLXXV   REM 
MCDLXXX   REM 
MCDLXXXV  REM 
MCDXC     REM 
MCDXCV    REM 
STOP
/bin/qbasic hack.bas
