#!/bin/bash

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Steven E. Bopp, Materials Science & Engineering-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-May 10, 2019:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-University of California, San Diego-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Begin Environment Directions-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Make Certain That All Pseudopotentials are Readable-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-and That all Executables are Called Below-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Make Certain That All Static Directory Links are Correct:-:-:-:-:-:-
#-:-:-:-:-:-:-:Install GNUPlot and Inkscape to Have Full Functionality:-:-:-:-:-:-
#-:-:-:-:-:-Make sure that the script is run as Bash and not anything else:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This code shows how to use pp.x to calculate the charge density of of TiN"

# set the needed environment variables
. /home/steven/Documents/qe-6.4.1/environment_variables

# required executables and pseudopotentials
BIN_LIST="pw.x pp.x plotrho.x bands.x plotband.x dos.x projwfc.x fs.x epsilon.x sumpdos.x"
PSEUDO_LIST="N_ONCV_PBE-1.0.upf Zr_ONCV_PBE-1.0.upf 
Zr.rel-pbesol-n-nc.UPF N.rel-pbesol-nc.UPF 
Zr.pz-spn-kjpaw_psl.1.0.0.UPF N.pz-n-kjpaw_psl.0.1.UPF 
Zr.pbesol-spn-kjpaw_psl.1.0.0.UPF N.pbesol-n-kjpaw_psl.1.0.0.UPF
Zr.pbesol-n-nc.UPF N.pbesol-nc.UPF
Ti.pbe-spn-kjpaw_psl.1.0.0.UPF N.pbe-n-kjpaw_psl.1.0.0.UPF
Ti.pz-n-nc.UPF N.pz-nc.UPF Zr.pz-n-nc.UPF
"

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
#INK_COMMAND=`which inkscape 2>/dev/null`
#if [ "$INK_COMMAND" = "" ]; then
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
# I usually just list all the .x files that I commonly use here so I don't have to worry about it in the future
       FS_COMMAND="$BIN_DIR/fs.x "
       PW_COMMAND="$PARA_PREFIX $BIN_DIR/pw.x $PARA_POSTFIX"
       PP_COMMAND="$PARA_PREFIX $BIN_DIR/pp.x $PARA_POSTFIX"
      DOS_COMMAND="$PARA_PREFIX $BIN_DIR/dos.x $PARA_POSTFIX"
      EPS_COMMAND="$PARA_PREFIX $BIN_DIR/epsilon.x $PARA_POSTFIX"
    BANDS_COMMAND="$PARA_PREFIX $BIN_DIR/bands.x $PARA_POSTFIX"
  PROJWFC_COMMAND="$PARA_PREFIX $BIN_DIR/projwfc.x $PARA_POSTFIX"
  PLOTRHO_COMMAND="$BIN_DIR/plotrho.x"
  SUMPDOS_COMMAND="$BIN_DIR/sumpdos.x"
 PLOTBAND_COMMAND="$BIN_DIR/plotband.x"

$ECHO
$ECHO "  running pw.x as:       $PW_COMMAND"
$ECHO "  running pp.x as:       $PP_COMMAND"
$ECHO "  running fs.x as:       $FS_COMMAND"
$ECHO "  running dos.x as:      $DOS_COMMAND"
$ECHO "  running bands.x as:    $BANDS_COMMAND"
$ECHO "  running gnuplot as:    $GP_COMMAND"
$ECHO "  running epsilon.x as:  $EPS_COMMAND"
$ECHO "  running plotrho.x as:  $PLOTRHO_COMMAND"
$ECHO "  running projwfc.x as:  $PROJWFC_COMMAND"
$ECHO "  running sumpdos.x as:  $SUMPDOS_COMMAND"
$ECHO "  running plotband.x as: $PLOTBAND_COMMAND"
$ECHO

START_TIME=$SECONDS # Begin elapsed time measurement (thanks Tom Anderson from StackOverflow)

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Begin Calculation Directions-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:Begin Self-Consistent Calculation:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-Calculating in order Gamma, X, W, K, L, Gamma-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:Use Norm-Conserving Pseudopotentials for epsilon.x to work-:-:-:-:-:-
#-:-:-:.in File partially Generated with Cif2Cell for the TiN B1 System-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-Gaussian Smearing is Suited for Metals:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-mv Smearing is the 1999 Marzari-Vanderbilt Cold Smearing:-:-:-:-:-:-:-
#-:-:wf_collect=.TRUE Sets the Number of Processors as Constant (for epsilin.x)-:-
#-:-:-:-:-:Make sure that ecutwfc and ecutrho are set to the same values:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:for the SCF and the NSCF calculations:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:Calculations done for an 8-atom cell with parameter half that of TiN-:-:-:-
#-:-:-:-:-:-:-:Since epsilon.x does not support the irreducible cell:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.scf.in << EOF
 &control
                 calculation = 'scf'
                restart_mode = 'from_scratch',
                      prefix = 'TiN'
                  pseudo_dir = '$PSEUDO_DIR/',
                      outdir = '$TMP_DIR/'
 /
 &system
                       ibrav = 0          ,
                   celldm(1) = 8.00299    , 
                         nat = 2          , 
                        ntyp = 2          ,
                     ecutwfc = 18.0       ,
                 occupations = 'smearing' ,
		    smearing = 'mv'       ,
                     degauss = 0.02       ,
 /
 &electrons
                    conv_thr =  1.0d-8
                 mixing_beta = 0.7
 /
 CELL_PARAMETERS {alat}
  0.500000000000000   0.500000000000000   0.000000000000000 
  0.500000000000000   0.000000000000000   0.500000000000000 
  0.000000000000000   0.500000000000000   0.500000000000000 
ATOMIC_SPECIES
   Ti   47.86700  Ti.pz-n-nc.UPF
   N    14.00650  N.pz-nc.UPF
 ATOMIC_POSITIONS {crystal}
   Ti   0.000000000000000   0.000000000000000   0.000000000000000 
   N    0.500000000000000   0.500000000000000   0.500000000000000
K_POINTS {automatic}
 20 20 20 0 0 0
EOF
$ECHO "  running an scf calculation for TiN...\c"
$PW_COMMAND < TiN.scf.in > TiN.scf.out
check_failure $?
$ECHO " done"

ELAPSED_TIME1=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME1 seconds"
 

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-Begin Charge Density Calculation with PP.x:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-Make sure that vectors e1 and e2 are orthogonal-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Make sure that isovalue is large enough for calculation-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:Set x0(...) at the origin of your crystal system-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Set iflag=2 for a 2D plot:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:Set output_format = 2 for plotrho:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:Set output_format = 3 for xcrysden-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-nfile = number of data files that you want to read:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Set plot_num=0 for electron (pseudo-)charge density-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.pp_rho.in << EOF
 &inputpp
                     prefix  = 'TiN'
                      outdir = '$TMP_DIR/'
                     filplot = 'TiNcharge'
                    plot_num = 0
 /
 &plot
                       nfile = 1
                   filepp(1) = 'TiNcharge'
                   weight(1) = 1.0
                       iflag = 2
               output_format = 2
                     fileout = 'TiN.rho.dat'
    e1(1) = 1.0, e1(2) = 1.0, e1(3) = 0.0,
    e2(1) = 0.0, e2(2) = 0.0, e2(3) = 1.0,
    nx=56, ny=40
 /
EOF
$ECHO "  running pp.x to make a 2-d plot of the charge density...\c"
$PP_COMMAND < TiN.pp_rho.in > TiN.pp_rho.out
check_failure $?
$ECHO " done"

ELAPSED_TIME2=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME2 seconds"


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Begin Plotrho Commands for electron (pseudo-)charge density-:-:-:-:-:-
#-:-:-First argument is the input file defined by fileout=... in the PP code:-:-:-
#-:-:-:-:-:-:-Second argument is the output file name that you want-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:Third argument gives y or n for logarithmic scale:-:-:-:-:-:-:-:-
#-:-Fourth argument gives rho_min rho_max and the number of coutour plot levels-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.plotrho.in << EOF
TiN.rho.dat
TiN.rho.ps
n
0 0.09 6
EOF

$ECHO "  running plotrho.x to generate rho.ps from the pp.x calculation...\c"
$PLOTRHO_COMMAND < TiN.plotrho.in > TiN.plotrho.out
$ECHO " done"

ELAPSED_TIME3=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME3 seconds"


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:Begin Charge Density Calculation with PP.x for (110)-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-Make sure that vectors e1 and e2 are orthogonal-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Make sure that isovalue is large enough for calculation-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:Set x0(...) at the origin of your crystal system-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Set iflag=2 for a 2D plot:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:Set output_format = 7 for gnuplot:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Set plot_num=0 for electron (pseudo-)charge density-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Increase nx and ny for a more smooth plot of TiN.charge.png-:-:-:-:-:-
#-:-:-:Original values for nx and ny are nx=141 and ny=100, coarse results:-:-:-:-
#-:-:-Changing the nx and ny may caouse the plor TiN.contour.ps to misbehave:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.pp_rho_110.in << EOF
 &inputpp
                     prefix  = 'TiN'
                      outdir = '$TMP_DIR/'
                     filplot = 'TiNcharge110'
                    plot_num = 0
 /
 &plot
                       nfile = 1
                   filepp(1) = 'TiNcharge110'
                   weight(1) = 1.0
                       iflag = 2
               output_format = 7
                     fileout = 'TiN.rho_110.dat'
    e1(1) = 1.0, e1(2) = 1.0, e1(3) = 0.0,
    e2(1) = 0.0, e2(2) = 0.0, e2(3) = 1.0,
    nx=56, ny=40
 /
EOF

$ECHO
$ECHO "  running pp.x for a (110) plot of the charge density...\c"
$PP_COMMAND < TiN.pp_rho_110.in > TiN.pp_rho_110.out
check_failure $?
$ECHO " done"

ELAPSED_TIME4=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME4 seconds"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Begin Plotting with GNUPlot for the (110) direction-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-set alat = celldm(1) for cubic systems, more complicated for non-cubic:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-set pm3d for 3D gnuplot-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:set xra and yra based on the direction that you are simulating-:-:-:-:-
#-:-for example in the (110) direction set x-range = [0:alat*sqrt(2)] for cubic-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-or set xra = [0:alat] for (100)-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-Do not forget to end the if statement with fi when editing this portion-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

if [ "$GP_COMMAND" = "" ]; then
    break
else
cat > gnuplot.tmp <<EOF
#!$GP_COMMAND
#
set term png font ",18" enh size 1000,707 
set pm3d
set palette model HSV functions gray*0.75, 1, 0.9
set view 0,0
#
alat=8.00299
set xra[0:1.4142136*alat]
set yra [0.:alat]
set size ratio 1./1.4142136
set xtics out nomirror
set ytics axis in offset -4.0,0 nomirror 
set label "r (a.u)" at 6.8,-2.2 center
set label "r (a.u)" at -1.7,5.0 rotate by 90 center
unset ztics
unset key
set colorbox 
#
set out 'TiN.charge110.png'
set title "TiN charge along 110"
splot 'TiN.rho_110.dat' u 1:2:3 w pm3d
EOF

$ECHO "  generating TiN.charge110.png...\c"
$GP_COMMAND < gnuplot.tmp
$ECHO " done"
#rm gnuplot.tmp

fi

ELAPSED_TIME5=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME5 seconds"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:Begin Charge Density Calculation with PP.x for (100)-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-Make sure that vectors e1 and e2 are orthogonal-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Make sure that isovalue is large enough for calculation-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:Set x0(...) at the origin of your crystal system-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Set iflag=2 for a 2D plot:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:Set output_format = 7 for gnuplot:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Set plot_num=0 for electron (pseudo-)charge density-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Increase nx and ny for a more smooth plot of TiN.charge.png-:-:-:-:-:-
#-:-:-:Original values for nx and ny are nx=141 and ny=100, coarse results:-:-:-:-
#-:-:-Changing the nx and ny may caouse the plor TiN.contour.ps to misbehave:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.pp_rho_100.in << EOF
 &inputpp
                     prefix  = 'TiN'
                      outdir = '$TMP_DIR/'
                     filplot = 'TiNcharge100'
                    plot_num = 0
 /
 &plot
                       nfile = 1
                   filepp(1) = 'TiNcharge100'
                   weight(1) = 1.0
                       iflag = 2
               output_format = 7
                     fileout = 'TiN.rho_100.dat'
    e1(1) = 1.0, e1(2) = 0.0, e1(3) = 0.0,
    e2(1) = 0.0, e2(2) = 0.0, e2(3) = 1.0,
    nx=56, ny=40
 /
EOF

$ECHO
$ECHO "  running pp.x for a (100) 2D plot of the charge density...\c"
$PP_COMMAND < TiN.pp_rho_100.in > TiN.pp_rho_100.out
check_failure $?
$ECHO " done"

ELAPSED_TIME6=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME6 seconds"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Begin Plotting with GNUPlot (Assuming That It is Installed)-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-set alat = celldm(1) for cubic systems, more complicated for non-cubic:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-set pm3d for 3D gnuplot-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:set xra and yra based on the direction that you are simulating-:-:-:-:-
#-:-for example in the (110) direction set x-range = [0:alat*sqrt(2)] for cubic-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-or set xra = [0:alat] for (100)-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-Do not forget to end the if statement with fi when editing this portion-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

if [ "$GP_COMMAND" = "" ]; then
    break
else
cat > gnuplot.tmp <<EOF
#!$GP_COMMAND
#
set term png font ",18" enh size 1000,1000 
set pm3d
set palette model HSV functions gray*0.75, 1, 0.9
set view 0,0
#
alat=8.00299
set xra[0:alat]
set yra [0.:alat]
set size ratio 1./1.
set xtics out nomirror
set ytics axis in offset -4.0,0 nomirror 
set label "r (a.u)" at 6.8,-2.2 center
set label "r (a.u)" at -1.7,5.0 rotate by 90 center
unset ztics
unset key
set colorbox 
#
set out 'TiN.charge100.png'
set title "TiN charge along 100"
splot 'TiN.rho_100.dat' u 1:2:3 w pm3d
EOF

$ECHO "  generating TiN.charge100.png...\c"
$GP_COMMAND < gnuplot.tmp
$ECHO " done"
#rm gnuplot.tmp

fi

ELAPSED_TIME7=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME7 seconds"


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:Begin Electron Localization Function Calculation with PP.x for (100)-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-Make sure that vectors e1 and e2 are orthogonal-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Make sure that isovalue is large enough for calculation-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:Set x0(...) at the origin of your crystal system-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Set iflag=2 for a 2D plot:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:Set output_format = 7 for gnuplot:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Set plot_num=0 for electron (pseudo-)charge density-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Increase nx and ny for a more smooth plot of TiN.charge.png-:-:-:-:-:-
#-:-:-:Original values for nx and ny are nx=141 and ny=100, coarse results:-:-:-:-
#-:-:-Changing the nx and ny may caouse the plor TiN.contour.ps to misbehave:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.pp_ELF_100.in << EOF
 &inputpp
                     prefix  = 'TiN'
                      outdir = '$TMP_DIR/'
                     filplot = 'TiN_ELF_100'
                    plot_num = 8
 /
 &plot
                       nfile = 1
                   filepp(1) = 'TiN_ELF_100'
                   weight(1) = 1.0
                       iflag = 2
               output_format = 7
                     fileout = 'TiN.ELF_100.dat'
    e1(1) = 1.0, e1(2) = 0.0, e1(3) = 0.0,
    e2(1) = 0.0, e2(2) = 0.0, e2(3) = 1.0,
    nx=56, ny=40
 /
EOF

$ECHO
$ECHO "  running pp.x  for a 2D plot of ELF...\c"
$PP_COMMAND < TiN.pp_ELF_100.in > TiN.pp_ELF_100.out
check_failure $?
$ECHO " done"

ELAPSED_TIME8=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME8 seconds"


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Begin Plotting with GNUPlot (ELF 100 2D Data Plot):-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-set alat = celldm(1) for cubic systems, more complicated for non-cubic:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-set pm3d for 3D gnuplot-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:set xra and yra based on the direction that you are simulating-:-:-:-:-
#-:-for example in the (110) direction set x-range = [0:alat*sqrt(2)] for cubic-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-or set xra = [0:alat] for (100)-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-Do not forget to end the if statement with fi when editing this portion-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

if [ "$GP_COMMAND" = "" ]; then
    break
else
cat > gnuplot.tmp <<EOF
#!$GP_COMMAND
#
set term png font ",18" enh size 1000,1000 
set pm3d
set palette model HSV functions gray*0.75, 1, 0.9
set view 0,0
#
alat=8.00299
set xra[0:alat]
set yra [0.:alat]
set size ratio 1./1.
set xtics out nomirror
set ytics axis in offset -4.0,0 nomirror 
set label "r (a.u)" at 6.8,-2.2 center
set label "r (a.u)" at -1.7,5.0 rotate by 90 center
unset ztics
unset key
set colorbox 
#
set out 'TiN.ELF100.png'
set title "TiN ELF along 100"
splot 'TiN.ELF_100.dat' u 1:2:3 w pm3d
EOF

$ECHO "  generating TiN.charge100.png...\c"
$GP_COMMAND < gnuplot.tmp
$ECHO " done"
#rm gnuplot.tmp

fi

ELAPSED_TIME9=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME9 seconds"


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:Begin Electron Pseudo Charge Density Calculation with PP.x for (100)-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-Make sure that vectors e1 and e2 are orthogonal-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Make sure that isovalue is large enough for calculation-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:Set x0(...) at the origin of your crystal system-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Set iflag=2 for a 2D plot:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:Set output_format = 7 for gnuplot:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Set plot_num=0 for electron (pseudo-)charge density-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Increase nx and ny for a more smooth plot of TiN.charge.png-:-:-:-:-:-
#-:-:-:Original values for nx and ny are nx=141 and ny=100, coarse results:-:-:-:-
#-:-:-Changing the nx and ny may caouse the plor TiN.contour.ps to misbehave:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > TiN.pp_EpCD_100.in << EOF
 &inputpp
                     prefix  = 'TiN'
                      outdir = '$TMP_DIR/'
                     filplot = 'TiN_EpCD_100'
                    plot_num = 0
	      spin_component = 0
 /
 &plot
                       nfile = 1
                   filepp(1) = 'TiN_EpCD_100'
                   weight(1) = 1.0
                       iflag = 2
               output_format = 7
                     fileout = 'TiN.EpCD_100.dat'
    e1(1) = 1.0, e1(2) = 0.0, e1(3) = 0.0,
    e2(1) = 0.0, e2(2) = 0.0, e2(3) = 1.0,
    nx=56, ny=40
 /
EOF

$ECHO
$ECHO "  running pp.x  for a 2D plot of ELF...\c"
$PP_COMMAND < TiN.pp_EpCD_100.in > TiN.pp_EpCD_100.out
check_failure $?
$ECHO " done"

ELAPSED_TIME10=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME10 seconds"


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Begin Plotting with GNUPlot (ELF 100 2D Data Plot):-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-set alat = celldm(1) for cubic systems, more complicated for non-cubic:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-set pm3d for 3D gnuplot-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:set xra and yra based on the direction that you are simulating-:-:-:-:-
#-:-for example in the (110) direction set x-range = [0:alat*sqrt(2)] for cubic-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-or set xra = [0:alat] for (100)-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-Do not forget to end the if statement with fi when editing this portion-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

if [ "$GP_COMMAND" = "" ]; then
    break
else
cat > gnuplot.tmp <<EOF
#!$GP_COMMAND
#
set term png font ",18" enh size 1000,1000 
set pm3d
set palette model HSV functions gray*0.75, 1, 0.9
set view 0,0
#
alat=8.00299
set xra[0:alat]
set yra [0.:alat]
set size ratio 1./1.
set xtics out nomirror
set ytics axis in offset -4.0,0 nomirror 
set label "r (a.u)" at 6.8,-2.2 center
set label "r (a.u)" at -1.7,5.0 rotate by 90 center
unset ztics
unset key
set colorbox 
#
set out 'TiN.EpCD100.png'
set title "TiN EpCD along 100"
splot 'TiN.EpCD_100.dat' u 1:2:3 w pm3d
EOF

$ECHO "  generating TiN.charge100.png...\c"
$GP_COMMAND < gnuplot.tmp
$ECHO " done"
#rm gnuplot.tmp

fi

ELAPSED_TIME11=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME11 seconds"


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-Begin Cleaning Temp Directory-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

$ECHO "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/TiN.*

$ECHO
$ECHO "$EXAMPLE_DIR: done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:End Elapsed Time Measurement-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME seconds"

echo "Job Completed"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:TiN_ELF_v1.sh:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Steven E. Bopp, Materials Science & Engineering-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-May 8, 2019-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-University of California, San Diego-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-End of File-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-















