rm trypass.bas
/bin/umodem trypass.bas STOP
X IF CHECKPASS(ARG(II), ARG(III)) THEN GOTO XXXX
XX PRINT "no" + CHR(X)
XXX END
XXXX PRINT "yes" + CHR(X)
STOP
/bin/qbasic trypass.bas
