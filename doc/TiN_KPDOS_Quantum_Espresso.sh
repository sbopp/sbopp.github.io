#!/bin/sh

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:TiN_KResolvedPDOS.sh-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Steven E. Bopp, Materials Science & Engineering-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:September 8, 2018:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-University of California, San Diego-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Begin Environment Directions-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Make Certain That All Pseudopotentials are Readable-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-and That all Executables are Called Below-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Make Certain That All Static Directory Links are Correct:-:-:-:-:-:-
#-:-:-:-:-:-:-:Install GNUPlot and Inkscape to Have Full Functionality:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-Fermi Energy for TiN is 17.2470 eV:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This code shows how to use pw.x and postprocessing codes to make a"
$ECHO "self consistent calculation for TiN, and to plot the band structure of"
$ECHO "TiN."

# set the needed environment variables
. /home/steven/Documents/Quantum_Espresso/qe-6.0/environment_variables

# required executables and pseudopotentials
BIN_LIST="pw.x pp.x plotrho.x bands.x plotband.x dos.x projwfc.x fs.x"
PSEUDO_LIST="Ti.pbe-spn-kjpaw_psl.1.0.0.UPF N.pbe-n-kjpaw_psl.1.0.0.UPF"
$ECHO
$ECHO "  executables directory: $BIN_DIR"
$ECHO "  pseudo directory:      $PSEUDO_DIR"
$ECHO "  temporary directory:   $TMP_DIR"
$ECHO
$ECHO "  checking that needed directories and files exist...\c"

# check for gnuplot
GP_COMMAND=`which gnuplot 2>/dev/null`
if [ "$GP_COMMAND" = "" ]; then
        $ECHO
        $ECHO "gnuplot not in PATH"
        $ECHO "Results will not be plotted"
fi

# check for inkscape
#Ink_COMMAND=`which inkscape 2>/dev/null`
#if [ "$Ink_COMMAND" = "" ]; then
#        $ECHO
#        $ECHO "inkscape not in PATH"
#        $ECHO "Results will not be converted to .png"
#fi

# check for directories
for DIR in "$BIN_DIR" "$PSEUDO_DIR" ; do
    if test ! -d $DIR ; then
        $ECHO
        $ECHO "ERROR: $DIR not existent or not a directory"
        $ECHO "Aborting"
        exit 1
    fi
done
for DIR in "$TMP_DIR" "$EXAMPLE_DIR/results" ; do
    if test ! -d $DIR ; then
        mkdir $DIR
    fi
done
cd $EXAMPLE_DIR/results

# check for executables
for FILE in $BIN_LIST ; do
    if test ! -x $BIN_DIR/$FILE ; then
        $ECHO
        $ECHO "ERROR: $BIN_DIR/$FILE not existent or not executable"
        $ECHO "Aborting"
        exit 1
    fi
done

# check for pseudopotentials
for FILE in $PSEUDO_LIST ; do
    if test ! -r $PSEUDO_DIR/$FILE ; then
       $ECHO
       $ECHO "Downloading $FILE to $PSEUDO_DIR...\c"
            $WGET $PSEUDO_DIR/$FILE $NETWORK_PSEUDO/$FILE 2> /dev/null
    fi
    if test $? != 0; then
        $ECHO
        $ECHO "ERROR: $PSEUDO_DIR/$FILE not existent or not readable"
        $ECHO "Aborting"
        exit 1
    fi
done
$ECHO " done"

# Executable locations, run instructions (stored as variables), and terminal text to be returned when run
       FS_COMMAND="$BIN_DIR/fs.x "
       PW_COMMAND="$PARA_PREFIX $BIN_DIR/pw.x $PARA_POSTFIX"
       PP_COMMAND="$PARA_PREFIX $BIN_DIR/pp.x $PARA_POSTFIX"
      DOS_COMMAND="$PARA_PREFIX $BIN_DIR/dos.x $PARA_POSTFIX"
    BANDS_COMMAND="$PARA_PREFIX $BIN_DIR/bands.x $PARA_POSTFIX"
  PROJWFC_COMMAND="$PARA_PREFIX $BIN_DIR/projwfc.x $PARA_POSTFIX"
  PLOTRHO_COMMAND="$BIN_DIR/plotrho.x"
 PLOTBAND_COMMAND="$BIN_DIR/plotband.x"
$ECHO
$ECHO "  running pw.x as:       $PW_COMMAND"
$ECHO "  running pp.x as:       $PP_COMMAND"
$ECHO "  running fs.x as:       $FS_COMMAND"
$ECHO "  running dos.x as:      $DOS_COMMAND"
$ECHO "  running bands.x as:    $BANDS_COMMAND"
$ECHO "  running gnuplot as:    $GP_COMMAND"
$ECHO "  running plotrho.x as:  $PLOTRHO_COMMAND"
$ECHO "  running projwfc.x as:  $PROJWFC_COMMAND"
$ECHO "  running plotband.x as: $PLOTBAND_COMMAND"
$ECHO

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Begin Calculation Directions-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:Begin Self-Consistent Calculation:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-Using kjpaw Pseudopotentials:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:.in File Generated with Cif2Cell for the TiN B1 System-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-Gaussian Smearing is Suited for Metals:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-mv Smearing is the 1999 Marzari-Vanderbilt Cold Smearing:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.scf.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='TiN'
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &SYSTEM
                       ibrav = 0          ,
                   celldm(1) = 8.00299    ,
                         nat = 2          ,
                        ntyp = 2          ,
                     ecutwfc = 40         ,
                     ecutrho = 160        ,
                 occupations = 'smearing' ,
		    smearing = 'mv'       ,
                     degauss = 0.02       ,
                       nspin = 2          ,
   starting_magnetization(1) = 1          ,
   starting_magnetization(2) = 1          ,
           tot_magnetization = -1         ,
 /
 &ELECTRONS
                    conv_thr = 1.0d-8     ,
                 mixing_beta = 0.7        ,
 /
 CELL_PARAMETERS {alat}
  0.500000000000000   0.500000000000000   0.000000000000000 
  0.500000000000000   0.000000000000000   0.500000000000000 
  0.000000000000000   0.500000000000000   0.500000000000000 
 ATOMIC_SPECIES
   Ti   47.86700  Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
   N    14.00650  N.pbe-n-kjpaw_psl.1.0.0.UPF
 ATOMIC_POSITIONS {crystal}
   Ti   0.000000000000000   0.000000000000000   0.000000000000000 
   N    0.500000000000000   0.500000000000000   0.500000000000000 

# k-space resolution ~0.2/A.
K_POINTS automatic
12 12 12  0 0 0

EOF
$ECHO "running the scf calculation...\c"
$PW_COMMAND < TiN.scf.in > TiN.scf.out
check_failure $?
$ECHO "done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:Begin Band Structure Calculation-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Using kjpaw Pseudopotentials-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:Calculating Along Delta and Sigma Lines:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.band.in << EOF
 &control
    calculation='bands'
    restart_mode='from_scratch',
    prefix='TiN',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &SYSTEM
                       ibrav = 2           ,
                   celldm(1) = 8.00299     ,
                         nat = 2           ,
                        ntyp = 2           ,
                     ecutwfc = 40          ,
                     ecutrho = 160         ,
                 occupations = 'smearing'  ,
		    smearing = 'mv'        ,
                     degauss = 0.02        ,
                       nspin = 2           ,
   starting_magnetization(1) = 1           ,
   starting_magnetization(2) = 1           ,
           tot_magnetization = -1          ,
 /
 &ELECTRONS
                    conv_thr = 1.0d-8      ,
                 mixing_beta = 0.7         ,
 /
 ATOMIC_SPECIES
   Ti   47.86700  Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
   N    14.00650  N.pbe-n-kjpaw_psl.1.0.0.UPF
 ATOMIC_POSITIONS {crystal}
   Ti   0.000000000000000   0.000000000000000   0.000000000000000 
   N    0.500000000000000   0.500000000000000   0.500000000000000 
K_POINTS
97
1.000000000 0.000000000 0.000000000 1
0.975000000 0.000000000 0.000000000 2
0.950000000 0.000000000 0.000000000 3
0.925000000 0.000000000 0.000000000 4
0.900000000 0.000000000 0.000000000 5
0.875000000 0.000000000 0.000000000 6
0.850000000 0.000000000 0.000000000 7
0.825000000 0.000000000 0.000000000 8
0.800000000 0.000000000 0.000000000 9
0.775000000 0.000000000 0.000000000 10
0.750000000 0.000000000 0.000000000 11
0.725000000 0.000000000 0.000000000 12
0.700000000 0.000000000 0.000000000 13
0.675000000 0.000000000 0.000000000 14
0.650000000 0.000000000 0.000000000 15
0.625000000 0.000000000 0.000000000 16
0.600000000 0.000000000 0.000000000 17
0.575000000 0.000000000 0.000000000 18
0.550000000 0.000000000 0.000000000 19
0.525000000 0.000000000 0.000000000 20
0.500000000 0.000000000 0.000000000 21
0.475000000 0.000000000 0.000000000 22
0.450000000 0.000000000 0.000000000 23
0.425000000 0.000000000 0.000000000 24
0.400000000 0.000000000 0.000000000 25
0.375000000 0.000000000 0.000000000 26
0.350000000 0.000000000 0.000000000 27
0.325000000 0.000000000 0.000000000 28
0.300000000 0.000000000 0.000000000 29
0.275000000 0.000000000 0.000000000 30
0.250000000 0.000000000 0.000000000 31
0.225000000 0.000000000 0.000000000 32
0.200000000 0.000000000 0.000000000 33
0.175000000 0.000000000 0.000000000 34
0.150000000 0.000000000 0.000000000 35
0.125000000 0.000000000 0.000000000 36
0.100000000 0.000000000 0.000000000 37
0.075000000 0.000000000 0.000000000 38
0.050000000 0.000000000 0.000000000 39
0.025000000 0.000000000 0.000000000 40
0.000000000 0.000000000 0.000000000 41
0.017857142 0.017857142 0.000000000 42
0.035714285 0.035714285 0.000000000 43
0.053571428 0.053571428 0.000000000 44
0.071428571 0.071428571 0.000000000 45
0.089285714 0.089285714 0.000000000 46
0.107142857 0.107142857 0.000000000 47
0.125000000 0.125000000 0.000000000 48
0.142857142 0.142857142 0.000000000 49
0.160714285 0.160714285 0.000000000 50
0.178571428 0.178571428 0.000000000 51
0.196428571 0.196428571 0.000000000 52
0.214285714 0.214285714 0.000000000 53
0.232142857 0.232142857 0.000000000 54
0.250000000 0.250000000 0.000000000 55
0.267857142 0.267857142 0.000000000 56
0.285714285 0.285714285 0.000000000 57
0.303571428 0.303571428 0.000000000 58
0.321428571 0.321428571 0.000000000 59
0.339285714 0.339285714 0.000000000 60
0.357142857 0.357142857 0.000000000 61
0.375000000 0.375000000 0.000000000 62
0.392857142 0.392857142 0.000000000 63
0.410714285 0.410714285 0.000000000 64
0.428571428 0.428571428 0.000000000 65
0.446428571 0.446428571 0.000000000 66
0.464285714 0.464285714 0.000000000 67
0.482142857 0.482142857 0.000000000 68
0.500000000 0.500000000 0.000000000 69
0.517857142 0.517857142 0.000000000 70
0.535714285 0.535714285 0.000000000 71
0.553571428 0.553571428 0.000000000 72
0.571428571 0.571428571 0.000000000 73
0.589285714 0.589285714 0.000000000 74
0.607142857 0.607142857 0.000000000 75
0.625000000 0.625000000 0.000000000 76
0.642857142 0.642857142 0.000000000 77
0.660714285 0.660714285 0.000000000 78
0.678571428 0.678571428 0.000000000 79
0.696428571 0.696428571 0.000000000 80
0.714285714 0.714285714 0.000000000 81
0.732142857 0.732142857 0.000000000 82
0.750000000 0.750000000 0.000000000 83
0.767857142 0.767857142 0.000000000 84
0.785714285 0.785714285 0.000000000 85
0.803571428 0.803571428 0.000000000 86
0.821428571 0.821428571 0.000000000 87
0.839285714 0.839285714 0.000000000 88
0.857142857 0.857142857 0.000000000 89
0.875000000 0.875000000 0.000000000 90
0.892857142 0.892857142 0.000000000 91
0.910714285 0.910714285 0.000000000 92
0.928571428 0.928571428 0.000000000 93
0.946428571 0.946428571 0.000000000 94
0.964285714 0.964285714 0.000000000 95
0.982142857 0.982142857 0.000000000 96
1.000000000 1.000000000 0.000000000 97

EOF
$ECHO "  running the band-structure calculation for TiN...\c"
$PW_COMMAND < TiN.band.in > TiN.band.out
check_failure $?
$ECHO " done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Begin K-resolved PDOS calculation along lines computed above:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:ngauss Is Gaussian Broadening Type, 0 Simple Default-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-0degauss Is Gaussian Broadening in Ry (not eV):-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:DeltaE Is Energy Grid Step Size in eV:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# 
cat > TiN.kpdos.in << EOF
 &projwfc
    outdir='$TMP_DIR/'
    prefix='TiN'
    ngauss=0, degauss=0.036748
    DeltaE=0.01
    kresolveddos=.true.
    filpdos='TiN.k'
 /
EOF
$ECHO "  running k-resolved PDOS calculation for TiN...\c"
$PROJWFC_COMMAND < TiN.kpdos.in > TiN.kpdos.out
check_failure $?
$ECHO " done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Begin Plotting with GNUPlot (Assuming That It is Installed)-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Begin Plotting with GNUPlot (Assuming That It is Installed)-:-:-:-:-:-
#-:-:-:-:-:--:-:-:-Don't forget to set the fermi energy (ef=...)-:-:-:-:-:--:-:-:-
#-:-:-:-Configured to Make Six Plots Since Multiplot Doesn't Like Many Plots:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
if [ "$GP_COMMAND" = "" ]; then
    break
else
cat > gnuplot.tmp <<EOF
#!$GP_COMMAND
#
set term png enh size 1000,500
set pm3d
set view 0,0
#
f(z)=z**(0.7)  # tune image contrast
ef=17.2470
#
unset xtics
set xtics out nomirror ("X" 1,"Gamma" 41,"K" 83, "X" 97)
set xra[1:97]
set label 1 "E-E_F(eV)" at 98,2.5
set ytics out nomirror
set yra [-20:7]
unset ztics
unset key
unset colorbox

#
set out 'kpdos_up1.png'
set origin 0,0
set size 1,1
set multiplot
dx=.1 ; dy=.30   # reduce margins
set title offset 0,-7
set size 1./3+1.4*dx,1.+2*dy
set origin 0./3-dx,0-dy
set title "Total DOS"
splot 'TiN.k.pdos_tot'		      u 1:(\$2-ef):(f(\$3)) w pm3d
set origin 1./3-dx,0-dy
set title "Ti_s1-DOS"
splot 'TiN.k.pdos_atm#1(Ti)_wfc#1(s)' u 1:(\$2-ef):(f(\$3)) w pm3d
set origin 2./3-dx,0-dy
set title "Ti_s2-DOS"
splot 'TiN.k.pdos_atm#1(Ti)_wfc#2(s)' u 1:(\$2-ef):(f(\$3)) w pm3d
unset multiplot

#
set out 'kpdos_up2.png'
set origin 0,0
set size 1,1
set multiplot
dx=.1 ; dy=.30   # reduce margins
set title offset 0,-7
set size 1./3+1.4*dx,1.+2*dy
set origin 0./3-dx,0-dy
set title " Total DOS"
splot 'TiN.k.pdos_tot'		      u 1:(\$2-ef):(f(\$3)) w pm3d
set origin 1./3-dx,0-dy
set title "Ti_p1-DOS"
splot 'TiN.k.pdos_atm#1(Ti)_wfc#3(p)' u 1:(\$2-ef):(f(\$3)) w pm3d
set origin 2./3-dx,0-dy
set title "Ti_d1-DOS"
splot 'TiN.k.pdos_atm#1(Ti)_wfc#4(d)' u 1:(\$2-ef):(f(\$3)) w pm3d
unset multiplot

#
set out 'kpdos_up3.png'
set origin 0,0
set size 1,1
set multiplot
dx=.1 ; dy=.30   # reduce margins
set title offset 0,-7
set size 1./3+1.4*dx,1.+2*dy
set origin 0./3-dx,0-dy
set title "Total DOS"
splot 'TiN.k.pdos_tot'		      u 1:(\$2-ef):(f(\$3)) w pm3d
set origin 1./3-dx,0-dy
set title "N_s1-DOS"
splot 'TiN.k.pdos_atm#2(N)_wfc#1(s)'  u 1:(\$2-ef):(f(\$3)) w pm3d
set origin 2./3-dx,0-dy
set title "N_p1-DOS"
splot 'TiN.k.pdos_atm#2(N)_wfc#2(p)'  u 1:(\$2-ef):(f(\$3)) w pm3d
unset multiplot

#
set out 'kpdos_dw1.png'
set origin 0,0
set size 1,1
set multiplot
dx=.1 ; dy=.30   # reduce margins
set title offset 0,-7
set size 1./3+1.4*dx,1.+2*dy
set origin 0./3-dx,0-dy
set title "Total DOS"
splot 'TiN.k.pdos_tot'		      u 1:(\$2-ef):(f(\$4)) w pm3d
set origin 1./3-dx,0-dy
set title "Ti_s1-DOS"
splot 'TiN.k.pdos_atm#1(Ti)_wfc#1(s)' u 1:(\$2-ef):(f(\$4)) w pm3d
set origin 2./3-dx,0-dy
set title "Ti_s2-DOS"
splot 'TiN.k.pdos_atm#1(Ti)_wfc#2(s)' u 1:(\$2-ef):(f(\$4)) w pm3d
unset multiplot

#
set out 'kpdos_dw2.png'
set origin 0,0
set size 1,1
set multiplot
dx=.1 ; dy=.30   # reduce margins
set title offset 0,-7
set size 1./3+1.4*dx,1.+2*dy
set origin 0./3-dx,0-dy
set title " Total DOS"
splot 'TiN.k.pdos_tot'		      u 1:(\$2-ef):(f(\$4)) w pm3d
set origin 1./3-dx,0-dy
set title "Ti_p1-DOS"
splot 'TiN.k.pdos_atm#1(Ti)_wfc#3(p)' u 1:(\$2-ef):(f(\$4)) w pm3d
set origin 2./3-dx,0-dy
set title "Ti_d1-DOS"
splot 'TiN.k.pdos_atm#1(Ti)_wfc#4(d)' u 1:(\$2-ef):(f(\$4)) w pm3d
unset multiplot

#
set out 'kpdos_dw3.png'
set origin 0,0
set size 1,1
set multiplot
dx=.1 ; dy=.30   # reduce margins
set title offset 0,-7
set size 1./3+1.4*dx,1.+2*dy
set origin 0./3-dx,0-dy
set title "Total DOS"
splot 'TiN.k.pdos_tot'		      u 1:(\$2-ef):(f(\$4)) w pm3d
set origin 1./3-dx,0-dy
set title "N_s1-DOS"
splot 'TiN.k.pdos_atm#2(N)_wfc#1(s)'  u 1:(\$2-ef):(f(\$4)) w pm3d
set origin 2./3-dx,0-dy
set title "N_p1-DOS"
splot 'TiN.k.pdos_atm#2(N)_wfc#2(p)'  u 1:(\$2-ef):(f(\$4)) w pm3d
unset multiplot

EOF
$ECHO
$ECHO "  plotting k-resolved DOS ...\c"
$GP_COMMAND < gnuplot.tmp
$ECHO " done"
rm gnuplot.tmp
fi

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-Begin Density of States Calculation-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:If One of the .out Files Returns an Error with 'Bad Fermi Energy':-:-:-:-
#-:-:-:-:-:-:-:-:-:-:Then Try Increasing nbnd to a larger number:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.dos.in << EOF
 &control
    calculation='nscf'
    prefix='TiN',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &SYSTEM
                       ibrav = 2            ,
                   celldm(1) = 8.00299      ,
                         nat = 2            ,
                        ntyp = 2            ,
                     ecutwfc = 40           ,
                     ecutrho = 160          ,
                 occupations = 'tetrahedra' ,
                     degauss = 0.02         ,
                       nspin = 2            ,
   starting_magnetization(1) = 1            ,
   starting_magnetization(2) = 1            ,
           tot_magnetization = -1           ,
                        nbnd = 10           ,
 /
 &electrons
                    conv_thr = 1.0d-8       ,
                 mixing_beta = 0.7          ,
 /
 ATOMIC_SPECIES
   Ti   47.86700  Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
   N    14.00650  N.pbe-n-kjpaw_psl.1.0.0.UPF
 ATOMIC_POSITIONS {crystal}
   Ti   0.000000000000000   0.000000000000000   0.000000000000000 
   N    0.500000000000000   0.500000000000000   0.500000000000000 
K_POINTS {automatic}
 12 12 12 0 0 0
EOF

cat > TiN.dos2.in << EOF
 &dos
    outdir='$TMP_DIR/'
    prefix='TiN'
    fildos='TiN.dos',
    Emin=5.0, Emax=25.0, DeltaE=0.1
 /
EOF

$ECHO "  running DOS calculation for TiN...\c"
$PW_COMMAND < TiN.dos.in > TiN.dos.out
check_failure $?
$DOS_COMMAND < TiN.dos2.in > TiN.dos2.out
check_failure $?
$ECHO " done"

cat > TiN.pdos.in << EOF
 &projwfc
    outdir='$TMP_DIR/'
    prefix='TiN'
    Emin=5.0, Emax=25.0, DeltaE=0.1
    ngauss=1, degauss=0.02
 /
EOF
$ECHO "  running PDOS calculation for TiN...\c"
$PROJWFC_COMMAND < TiN.pdos.in > TiN.pdos.out
check_failure $?
$ECHO " done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:Update User for Fermi Surface Section:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
$ECHO
$ECHO "  Beginning Calculations for the Fermi Surface plot, "
$ECHO "  Spin-Polarized case..."

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:Begin Self-Consistent Calculation for the Spin-PolariTid (SP) Case-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.scf_SP.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='TiN',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &SYSTEM
                       ibrav = 2                  ,
                   celldm(1) = 8.00299            ,
                         nat = 2                  ,
                        ntyp = 2                  ,
                     ecutwfc = 40                 ,
                     ecutrho = 160                ,
                 occupations = 'smearing'         ,
                    smearing = 'methfessel-paxton',
                     degauss = 0.02               ,
                       nspin = 2                  ,
   starting_magnetization(1) = 1                  ,
   starting_magnetization(2) = 1                  ,
           tot_magnetization = -1                 ,
 /
 &electrons
                    conv_thr = 1.0d-8             ,
                 mixing_beta = 0.7                ,
 /
 ATOMIC_SPECIES
   Ti   47.86700  Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
   N    14.00650  N.pbe-n-kjpaw_psl.1.0.0.UPF
 ATOMIC_POSITIONS {crystal}
   Ti   0.000000000000000   0.000000000000000   0.000000000000000 
   N    0.500000000000000   0.500000000000000   0.500000000000000 
K_POINTS {automatic}
 8 8 8 0 0 0 
EOF
$ECHO "  running the scf calculation for the spin-polarized case ... \c"
$PW_COMMAND < TiN.scf_SP.in > TiN.scf0.SP.out
check_failure $?
$ECHO " done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-Begin Fermi Surface Calculation for the Spin-PolariTid (SP) Case:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-If This Section Halts Unexpectedly with The Following Errors:-:-:-:-:-
#-:-:-:-:1 Error condition encountered during test: exit status = 1 Aborting:-:-:-
#-:-:-:-:-:-:-:-:-:-:Then Try Increasing nbnd to a larger number:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.fs_SP.in << EOF
 &control
    calculation='nscf'
    prefix='TiN',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &SYSTEM
                       ibrav = 2            ,
                   celldm(1) = 8.00299      ,
                         nat = 2            ,
                        ntyp = 2            ,
                     ecutwfc = 40           ,
                     ecutrho = 160          ,
                 occupations = 'tetrahedra' ,
                        nbnd = 10           ,
                       nspin = 2            ,
   starting_magnetization(1) = 1            ,
   starting_magnetization(2) = 1            ,
           tot_magnetization = -1           ,
 /
 &electrons
                    conv_thr = 1.0d-8       ,
                 mixing_beta = 0.7          ,
 /
 ATOMIC_SPECIES
   Ti   47.86700  Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
   N    14.00650  N.pbe-n-kjpaw_psl.1.0.0.UPF
 ATOMIC_POSITIONS {crystal}
   Ti   0.000000000000000   0.000000000000000   0.000000000000000 
   N    0.500000000000000   0.500000000000000   0.500000000000000 
K_POINTS {automatic}
16 16 16 0 0 0 
EOF

$ECHO "  running the Fermi Surface calculation ... \c"
$PW_COMMAND   < TiN.fs_SP.in > TiN.fs_SP.out 
check_failure $?
$ECHO " done"

cat > FS.in <<EOF
&fermi
    outdir='$TMP_DIR/'
    prefix='TiN'
/
EOF
$FS_COMMAND < FS.in > FS.out
check_failure $?

$ECHO
$ECHO "  Use 'xcrysden --bxsf results/TiN_fsXX.bxsf', XX=up,dw to plot the Fermi Surface\c"
$ECHO " done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:Update User for Fermi Surface Section:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
$ECHO
$ECHO "  Beginning Calculations for the Fermi Surface plot, "
$ECHO "  Non-Spin-Polarized case..."

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-Begin Self-Consistent Calculation for the Non-Spin-PolariTid (NSP) Case-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.scf_NSP.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='TiN',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &SYSTEM
                       ibrav = 2                  ,
                   celldm(1) = 8.00299            ,
                         nat = 2                  ,
                        ntyp = 2                  ,
                     ecutwfc = 40                 ,
                     ecutrho = 160                ,
                 occupations = 'smearing'         ,
                    smearing = 'methfessel-paxton',
                     degauss = 0.02               ,
                        nbnd = 10                 ,
                       nspin = 2                  ,
   starting_magnetization(1) = 1                  ,
   starting_magnetization(2) = 1                  ,
           tot_magnetization = -1                 ,
 /
 &electrons
                    conv_thr = 1.0d-8             ,
                 mixing_beta = 0.7                ,
 /
 ATOMIC_SPECIES
   Ti   47.86700  Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
   N    14.00650  N.pbe-n-kjpaw_psl.1.0.0.UPF
 ATOMIC_POSITIONS {crystal}
   Ti   0.000000000000000   0.000000000000000   0.000000000000000 
   N    0.500000000000000   0.500000000000000   0.500000000000000 
K_POINTS {automatic}
 8 8 8 0 0 0 
EOF
$ECHO "  running the scf calculation  non spin-polarized case ... \c"
$PW_COMMAND < TiN.scf_NSP.in > TiN.scf0.NSP.out
check_failure $?
$ECHO " done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-Begin Fermi Surface Calculation for the Spin-PolariTid (NSP) Case-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-If This Section Halts Unexpectedly with The Following Errors:-:-:-:-:-
#-:-:-:-:1 Error condition encountered during test: exit status = 1 Aborting:-:-:-
#-:-:-:-:-:-:-:-:-:-:Then Try Increasing nbnd to a larger number:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.fs_NSP.in << EOF
 &control
    calculation='nscf'
    prefix='TiN',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &SYSTEM
                       ibrav = 2            ,
                   celldm(1) = 8.00299      ,
                         nat = 2            ,
                        ntyp = 2            ,
                     ecutwfc = 40           ,
                     ecutrho = 160          ,
                 occupations = 'tetrahedra' ,
                        nbnd = 10           ,
 /
 &electrons
                    conv_thr = 1.0d-8       ,
                 mixing_beta = 0.7          ,
 /
 ATOMIC_SPECIES
   Ti   47.86700  Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
   N    14.00650  N.pbe-n-kjpaw_psl.1.0.0.UPF
 ATOMIC_POSITIONS {crystal}
   Ti   0.000000000000000   0.000000000000000   0.000000000000000 
   N    0.500000000000000   0.500000000000000   0.500000000000000 
K_POINTS automatic
16 16 16 0 0 0 
EOF

$ECHO "  running the Fermi Surface calculation ... \c"
$PW_COMMAND   < TiN.fs_NSP.in > TiN.fs_NSP.out 
check_failure $?
$ECHO " done"

cat > FS.in <<EOF
&fermi
    outdir='$TMP_DIR/'
    prefix='TiN'
/
EOF
$FS_COMMAND < FS.in > FS.out
check_failure $?

$ECHO
$ECHO "  Use 'xcrysden --bxsf results/TiN_fs.bxsf' to plot the Fermi Surface\c"
$ECHO " done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-Begin Cleaning Temp Directory-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

$ECHO "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/TiN.*

$ECHO
$ECHO "$EXAMPLE_DIR: done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:TiN_KResolvedPDOS.sh-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Steven E. Bopp, Materials Science & Engineering-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:September 8, 2018:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-University of California, San Diego-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-End of File-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-










