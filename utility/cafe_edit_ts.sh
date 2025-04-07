#!/bin/bash
SFILE=1
GPLOT=1
SGRAPH=1

export UNIT_NAME="#all"
FNAME=$1_all

XFIELD=1
YFIELD=6
XAXIS="Time"
YAXIS="Energy"
XTICS=
YTICS=
XRANGE="[:]"
YRANGE="[:]"

if test $SFILE -eq 1
    then
    uni2=`awk '/#unit-unit/{print "2"}' $1`
    if test $uni2 -eq 2
	then
	awk '/^#all      /{print $0}' $1 | sed 's/#.........//' > $FNAME
    else
	awk '/^#all /{print $0}' $1 | sed 's/#.........//' > $FNAME
    fi
fi


if test $GPLOT -eq 1
    then
    echo "set xlabel \"$XAXIS\"; \
          set ylabel \"$YAXIS\"; \
          set xtics $XTICS; \
          set ytics $YTICS; \
          plot $XRANGE $YRANGE \"$FNAME\" u $XFIELD:$YFIELD w l" | gnuplot -persist
fi


if test $SGRAPH -eq 1
    then
#    colors=( red orange yellow green black blue purple )
    echo "set terminal postscript eps color enhanced; \
          unset key; set output 'test.eps'; \
          set xlabel \"$XAXIS\"; \
          set ylabel \"$YAXIS\"; \
          set xtics $XTICS; \
          set ytics $YTICS; \
          plot $XRANGE $YRANGE '$FNAME' u $XFIELD:$YFIELD w l lt 1" | gnuplot
#    plot '$1_all' w l lt 1 lc rgb '${colors[0]}'" >> tmp
fi
