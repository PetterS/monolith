#!/bin/sh
debug=no
log=undefined
check=undefined
static=no
usage () {
cat<<EOF
usage: configure.sh [<option> ...]

-h|--help    print this command line option summary
-g|--debug   include debugging code (disabled by default)
-c|--check   include checking code (default for '-g')
-l|--log     include logging code (default for '-g')
-s|--static  static compilation
EOF
}
while [ $# -gt 0 ]
do
  case $1 in
    -h|--hepp) usage; exit 0;;
    -g|--debug) debug=yes;;
    -c|--check) check=yes;;
    -l|--log) log=yes;;
    -s|-static|--static) static=yes;;
    *) echo "*** configure.sh: invalid option '$1' (try '-h')"; exit 1;;
  esac
  shift
done
[ $check = undefined ] && check=$debug
[ $log = undefined ] && log=$debug
CXX=g++
CXXFLAGS="-Wall"
LDFLAGS=""
if [ $debug = yes ]
then
  CXXFLAGS="$CXXFLAGS -g3"
else
  CXXFLAGS="$CXXFLAGS -O3"
fi
[ $static = yes ] && CXXFLAGS="$CXXFLAGS -static"
[ $log = no ] && CXXFLAGS="$CXXFLAGS -DNLOG"
if [ $check = yes ]
then
  if [ -d ../lingeling -a ../lingeling/liblgl.a -a ../lingeling/lglib.h ]
  then
    CXXFLAGS="$CXXFLAGS -I../lingeling"
    LDFLAGS="-L../lingeling -llgl"
  else
    CXXFLAGS="$CXXFLAGS -DNLGL"
  fi
else
  CXXFLAGS="$CXXFLAGS -DNLGL -DNDEBUG"
fi
echo "$CXX $CXXFLAGS $LDFLAGS"
rm -f makefile
sed \
  -e "s,@CXX@,$CXX," \
  -e "s,@CXXFLAGS@,$CXXFLAGS," \
  -e "s,@LDFLAGS@,$LDFLAGS," \
  makefile.in > makefile

