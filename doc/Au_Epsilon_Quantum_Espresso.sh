#!/bin/bash

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:Au_Epsilon.sh:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Steven E. Bopp, Materials Science & Engineering-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-September 29, 2018:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-University of California, San Diego-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:Begin Environment Directions-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Make Certain That All Pseudopotentials are Readable-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-and That all Executables are Called Below-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-Make Certain That All Static Directory Links are Correct:-:-:-:-:-:-
#-:-:-:-:-:-:-:Install GNUPlot and Inkscape to Have Full Functionality:-:-:-:-:-:-
#-:-:-:-:-:-:Make that the shell is set to BASH, remove timer otherwise-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This code shows how to use epsilon.x to calculate the permittivity of Au"
$ECHO "Permittivity is calculated after SCF and NSCF calculations"
$ECHO "and then plotted with literature values using Gnuplot"

# set the needed environment variables
. /home/steven/Documents/Quantum_Espresso/qe-6.0/environment_variables

# required executables and pseudopotentials
BIN_LIST="pw.x pp.x plotrho.x bands.x plotband.x dos.x projwfc.x fs.x epsilon.x"
PSEUDO_LIST="
Au.pbe-n-nc.UPF 
Au_ONCV_PBE-1.0.upf 
Au.rel-pz-n-nc.UPF 
Au.pz-rrkjus_aewfc.UPF 
Au.bp-n-nc.UPF 
Au.rel-pbesol-n-nc.UPF"

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
       FS_COMMAND="$BIN_DIR/fs.x "
       PW_COMMAND="$PARA_PREFIX $BIN_DIR/pw.x $PARA_POSTFIX"
       PP_COMMAND="$PARA_PREFIX $BIN_DIR/pp.x $PARA_POSTFIX"
      DOS_COMMAND="$PARA_PREFIX $BIN_DIR/dos.x $PARA_POSTFIX"
      EPS_COMMAND="$PARA_PREFIX $BIN_DIR/epsilon.x $PARA_POSTFIX"
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
$ECHO "  running epsilon.x as:  $EPS_COMMAND"
$ECHO "  running plotrho.x as:  $PLOTRHO_COMMAND"
$ECHO "  running projwfc.x as:  $PROJWFC_COMMAND"
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
#-:-:-:.in File partially Generated with Cif2Cell for the Au B1 System-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-Gaussian Smearing is Suited for Metals:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-mv Smearing is the 1999 Marzari-Vanderbilt Cold Smearing:-:-:-:-:-:-:-
#-:-:wf_collect=.TRUE Sets the Number of Processors as Constant (for epsilin.x)-:-
#-:-:-:-:-:Make sure that ecutwfc and ecutrho are set to the same values:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:for the SCF and the NSCF calculations:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:Calculations done for a single-atom cell:-:-:-:--:-:-:-:-:-:-
#-:-:-:-:-:-:-:Since epsilon.x does not support the irreducible cell:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > Au.scf.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='Au'
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    wf_collect=.TRUE,
 /
 &SYSTEM
                       ibrav = 0          ,
                   celldm(1) = 5.44951    ,
                         nat = 1          ,
                        ntyp = 1          ,
                     ecutwfc = 65         ,
                     ecutrho = 260        ,
                 occupations = 'smearing' ,
		    smearing = 'mv'       ,
                     degauss = 0.02       ,
 /
 &ELECTRONS
                    conv_thr = 1.0d-8     ,
                 mixing_beta = 0.7        ,
 /
CELL_PARAMETERS {alat}
  0.577350269189626   0.000000000000000   0.816496580927725 
 -0.288675134594812  -0.500000000000000   0.816496580927725 
 -0.288675134594812   0.500000000000000   0.816496580927725 
ATOMIC_SPECIES
  Au  196.96600  Au.rel-pbesol-n-nc.UPF 
ATOMIC_POSITIONS {crystal}
Au   0.000000000000000   0.000000000000000   0.000000000000000 

 K_POINTS crystal
          60
        0.0000000000    0.0000000000    0.0000000000    1.0
        0.0000000000    0.0277777778    0.0277777778    1.0
        0.0000000000    0.0555555556    0.0555555556    1.0
        0.0000000000    0.0833333333    0.0833333333    1.0
        0.0000000000    0.1111111111    0.1111111111    1.0
        0.0000000000    0.1388888889    0.1388888889    1.0
        0.0000000000    0.1666666667    0.1666666667    1.0
        0.0000000000    0.1944444444    0.1944444444    1.0
        0.0000000000    0.2222222222    0.2222222222    1.0
        0.0000000000    0.2500000000    0.2500000000    1.0
        0.0000000000    0.2777777778    0.2777777778    1.0
        0.0000000000    0.3055555556    0.3055555556    1.0
        0.0000000000    0.3333333333    0.3333333333    1.0
        0.0000000000    0.3611111111    0.3611111111    1.0
        0.0000000000    0.3888888889    0.3888888889    1.0
        0.0000000000    0.4166666667    0.4166666667    1.0
        0.0000000000    0.4444444444    0.4444444444    1.0
        0.0000000000    0.4722222222    0.4722222222    1.0
        0.0000000000    0.5000000000    0.5000000000    1.0
        0.0277777778    0.5277777778    0.5000000000    1.0
        0.0555555556    0.5555555556    0.5000000000    1.0
        0.0833333333    0.5833333333    0.5000000000    1.0
        0.1111111111    0.6111111111    0.5000000000    1.0
        0.1388888889    0.6388888889    0.5000000000    1.0
        0.1666666667    0.6666666667    0.5000000000    1.0
        0.1944444444    0.6944444444    0.5000000000    1.0
        0.2222222222    0.7222222222    0.5000000000    1.0
        0.2500000000    0.7500000000    0.5000000000    1.0
        0.2708333333    0.7500000000    0.4791666667    1.0
        0.2916666667    0.7500000000    0.4583333333    1.0
        0.3125000000    0.7500000000    0.4375000000    1.0
        0.3333333333    0.7500000000    0.4166666667    1.0
        0.3541666667    0.7500000000    0.3958333333    1.0
        0.3750000000    0.7500000000    0.3750000000    1.0
        0.3863636364    0.7272727273    0.3863636364    1.0
        0.3977272727    0.7045454545    0.3977272727    1.0
        0.4090909091    0.6818181818    0.4090909091    1.0
        0.4204545455    0.6590909091    0.4204545455    1.0
        0.4318181818    0.6363636364    0.4318181818    1.0
        0.4431818182    0.6136363636    0.4431818182    1.0
        0.4545454545    0.5909090909    0.4545454545    1.0
        0.4659090909    0.5681818182    0.4659090909    1.0
        0.4772727273    0.5454545455    0.4772727273    1.0
        0.4886363636    0.5227272727    0.4886363636    1.0
        0.5000000000    0.5000000000    0.5000000000    1.0
        0.4666666667    0.4666666667    0.4666666667    1.0
        0.4333333333    0.4333333333    0.4333333333    1.0
        0.4000000000    0.4000000000    0.4000000000    1.0
        0.3666666667    0.3666666667    0.3666666667    1.0
        0.3333333333    0.3333333333    0.3333333333    1.0
        0.3000000000    0.3000000000    0.3000000000    1.0
        0.2666666667    0.2666666667    0.2666666667    1.0
        0.2333333333    0.2333333333    0.2333333333    1.0
        0.2000000000    0.2000000000    0.2000000000    1.0
        0.1666666667    0.1666666667    0.1666666667    1.0
        0.1333333333    0.1333333333    0.1333333333    1.0
        0.1000000000    0.1000000000    0.1000000000    1.0
        0.0666666667    0.0666666667    0.0666666667    1.0
        0.0333333333    0.0333333333    0.0333333333    1.0
        0.0000000000    0.0000000000    0.0000000000    1.0

EOF
$ECHO "running the scf calculation...\c"
$PW_COMMAND < Au.scf.in > Au.scf.out
check_failure $?
$ECHO "done"

ELAPSED_TIME1=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME1 seconds"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:Begin Non-Self-Consistent Calculation:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:Use Norm-Conserving Pseudopotentials for epsilon.x to work-:-:-:-:-:-	
#-:-:-:-:-:-:-:-:-:Run NSCF to calculate unoccupied bands because-:-:-:-:-:-:-:-:-
#-:-:-:-:-metals may not have a gap in the range over which QE calculates-:-:-:-:-
#-:-:-:-:-:-Reduced number of K_POINTS to reduce computational intensity:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > Au.nscf.in << EOF
 &control
    calculation = 'nscf'
    restart_mode='from_scratch',
    prefix='Au',
    tprnfor = .true.
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &system
    ibrav=  2, celldm(1) =5.44951, nat=  1, ntyp= 1,
    lspinorb=.true.,
    noncolin=.true.,
    starting_magnetization=0.0,
    occupations='smearing',
    degauss=0.02,
    smearing='mp',
    ecutwfc =65,
    ecutrho =260,
 /
 &electrons
    mixing_beta = 0.7,
    conv_thr =  1.0d-8
 /
ATOMIC_SPECIES
  Au  196.96600  Au.rel-pbesol-n-nc.UPF 
ATOMIC_POSITIONS {crystal}
  Au   0.000000000000000   0.000000000000000   0.000000000000000 

 K_POINTS crystal
          60
        0.0000000000    0.0000000000    0.0000000000    1.0
        0.0000000000    0.0277777778    0.0277777778    1.0
        0.0000000000    0.0555555556    0.0555555556    1.0
        0.0000000000    0.0833333333    0.0833333333    1.0
        0.0000000000    0.1111111111    0.1111111111    1.0
        0.0000000000    0.1388888889    0.1388888889    1.0
        0.0000000000    0.1666666667    0.1666666667    1.0
        0.0000000000    0.1944444444    0.1944444444    1.0
        0.0000000000    0.2222222222    0.2222222222    1.0
        0.0000000000    0.2500000000    0.2500000000    1.0
        0.0000000000    0.2777777778    0.2777777778    1.0
        0.0000000000    0.3055555556    0.3055555556    1.0
        0.0000000000    0.3333333333    0.3333333333    1.0
        0.0000000000    0.3611111111    0.3611111111    1.0
        0.0000000000    0.3888888889    0.3888888889    1.0
        0.0000000000    0.4166666667    0.4166666667    1.0
        0.0000000000    0.4444444444    0.4444444444    1.0
        0.0000000000    0.4722222222    0.4722222222    1.0
        0.0000000000    0.5000000000    0.5000000000    1.0
        0.0277777778    0.5277777778    0.5000000000    1.0
        0.0555555556    0.5555555556    0.5000000000    1.0
        0.0833333333    0.5833333333    0.5000000000    1.0
        0.1111111111    0.6111111111    0.5000000000    1.0
        0.1388888889    0.6388888889    0.5000000000    1.0
        0.1666666667    0.6666666667    0.5000000000    1.0
        0.1944444444    0.6944444444    0.5000000000    1.0
        0.2222222222    0.7222222222    0.5000000000    1.0
        0.2500000000    0.7500000000    0.5000000000    1.0
        0.2708333333    0.7500000000    0.4791666667    1.0
        0.2916666667    0.7500000000    0.4583333333    1.0
        0.3125000000    0.7500000000    0.4375000000    1.0
        0.3333333333    0.7500000000    0.4166666667    1.0
        0.3541666667    0.7500000000    0.3958333333    1.0
        0.3750000000    0.7500000000    0.3750000000    1.0
        0.3863636364    0.7272727273    0.3863636364    1.0
        0.3977272727    0.7045454545    0.3977272727    1.0
        0.4090909091    0.6818181818    0.4090909091    1.0
        0.4204545455    0.6590909091    0.4204545455    1.0
        0.4318181818    0.6363636364    0.4318181818    1.0
        0.4431818182    0.6136363636    0.4431818182    1.0
        0.4545454545    0.5909090909    0.4545454545    1.0
        0.4659090909    0.5681818182    0.4659090909    1.0
        0.4772727273    0.5454545455    0.4772727273    1.0
        0.4886363636    0.5227272727    0.4886363636    1.0
        0.5000000000    0.5000000000    0.5000000000    1.0
        0.4666666667    0.4666666667    0.4666666667    1.0
        0.4333333333    0.4333333333    0.4333333333    1.0
        0.4000000000    0.4000000000    0.4000000000    1.0
        0.3666666667    0.3666666667    0.3666666667    1.0
        0.3333333333    0.3333333333    0.3333333333    1.0
        0.3000000000    0.3000000000    0.3000000000    1.0
        0.2666666667    0.2666666667    0.2666666667    1.0
        0.2333333333    0.2333333333    0.2333333333    1.0
        0.2000000000    0.2000000000    0.2000000000    1.0
        0.1666666667    0.1666666667    0.1666666667    1.0
        0.1333333333    0.1333333333    0.1333333333    1.0
        0.1000000000    0.1000000000    0.1000000000    1.0
        0.0666666667    0.0666666667    0.0666666667    1.0
        0.0333333333    0.0333333333    0.0333333333    1.0
        0.0000000000    0.0000000000    0.0000000000    1.0

EOF
$ECHO "  running the non-scf calculation for Au with norm-conserving PP...\c"
$PW_COMMAND < Au.nscf.in > Au.nscf.out
check_failure $?
$ECHO " done"

ELAPSED_TIME2=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME2 seconds"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-Begin epsilon.x Calculation-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:Use Norm-Conserving Pseudopotentials for epsilon.x to work-:-:-:-:-:-	
#-:-:-:-:-:-:-:-:-:-:Calculating from the SCF and NSCF aboove-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-intrasmear is the broadening parameter for the intraband:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-i.e. metal Drude like term (again in eV):-:-:-:-:-:-:-:-:-:-
#-:intersmear is the broadening parameter (in eV) for the interband contribution:-
#-:-:-:-:-:-:-:-:-:-:-:Intersmear default is 0.136d0 (eV)-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-Interband transition threshold ~2.5 for Au:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-nw is the number of points in the frequency mesh:-:-:-:-:-:-:-:-
#-:-:-:-:-Calculation performed on interval [-wmax,wmax], wmax units of eV:-:-:-:-
#-Shift is number of eV for rigid shift of the imaginary part of the wavefunction-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-If epsilon.x throws an error concerning division by zero:-:-:-:-:-:-
#-:-:-:then it is likely that intersmear is too small giving and infinities-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
cat > Au.eps.in << EOF
 &inputpp
    outdir='$TMP_DIR/'
    prefix='Au'
    calculation='eps'
 /
    &energy_grid
    smeartype='lorentz'
    intersmear=2.30d0
    intrasmear=0.00d0
    wmax=6.5d0
    wmin=0.0d0
    nw=500
    shift=0.0d0
 /
EOF
$ECHO "running the permittivity calculation...\c"
$EPS_COMMAND < Au.eps.in > Au.eps.out
check_failure $?
$ECHO "done"

ELAPSED_TIME3=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME3 seconds"


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:Create Au_Permittivity.dat-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:Values taken from:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:P. B. Johnson and R. W. Christy. Optical constants of the noble metals,:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-Phys. Rev. B 6, 4370-4379 (1972):-:-:-:-:-:-:-:-:-:-:-:-
#-:-:Values sorted as Wavelength (nm), n, k, energy (eV), real, and imaginary-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
$ECHO "Creating Au_Permittivity.dat from Johnson Christy Data...\c"

cat > Au.dat << EOF

1.93700000 0.92000000 13.78000000 0.64016520 -189.04200000 25.35520000
1.61000000 0.56000000 11.21000000 0.77018634 -125.35050000 12.55520000
1.39300000 0.43000000 9.51900000 0.89016511 -90.42646100 8.18634000
1.21600000 0.35000000 8.14500000 1.01973684 -66.21852500 5.70150000
1.08800000 0.27000000 7.15000000 1.13970588 -51.04960000 3.86100000
0.98400000 0.22000000 6.35000000 1.26016260 -40.27410000 2.79400000
0.89200000 0.17000000 5.66300000 1.39013453 -32.04066900 1.92542000
0.82110000 0.16000000 5.08300000 1.51016929 -25.81128900 1.62656000
0.75600000 0.14000000 4.54200000 1.64021164 -20.61016400 1.27176000
0.70450000 0.13000000 4.10300000 1.76011356 -16.81770900 1.06678000
0.65950000 0.14000000 3.69700000 1.88021228 -13.64820900 1.03516000
0.61680000 0.21000000 3.27200000 2.01037613 -10.66188400 1.37424000
0.58210000 0.29000000 2.86300000 2.13021818 -8.11266900 1.66054000
0.54860000 0.43000000 2.45500000 2.26029894 -5.84212500 2.11130000
0.52090000 0.62000000 2.08100000 2.38049530 -3.94616100 2.58044000
0.49590000 1.04000000 1.83300000 2.50050413 -2.27828900 3.81264000
0.47140000 1.31000000 1.84900000 2.63046245 -1.70270100 4.84438000
0.45090000 1.38000000 1.91400000 2.75005544 -1.75899600 5.28264000
0.43050000 1.45000000 1.94800000 2.88037166 -1.69220400 5.64920000
0.41330000 1.46000000 1.95800000 3.00024195 -1.70216400 5.71736000
0.39740000 1.47000000 1.95200000 3.12028183 -1.64940400 5.73888000
0.38150000 1.46000000 1.93300000 3.25032765 -1.60488900 5.64436000
0.36790000 1.48000000 1.89500000 3.37048111 -1.40062500 5.60920000
0.35420000 1.50000000 1.86600000 3.50084698 -1.23195600 5.59800000
0.34250000 1.48000000 1.87100000 3.62043796 -1.31024100 5.53816000
0.33150000 1.48000000 1.88300000 3.74057315 -1.35528900 5.57368000
0.32040000 1.54000000 1.89800000 3.87016230 -1.23080400 5.84584000
0.31070000 1.53000000 1.89300000 3.99098809 -1.24254900 5.79258000
0.30090000 1.53000000 1.88900000 4.12097042 -1.22742100 5.78034000
0.29240000 1.49000000 1.87800000 4.24076607 -1.30678400 5.59644000
0.28440000 1.47000000 1.86900000 4.36005626 -1.33226100 5.49486000
0.27610000 1.43000000 1.84700000 4.49112640 -1.36650900 5.28242000
0.26890000 1.38000000 1.80300000 4.61137970 -1.34640900 4.97628000
0.26160000 1.35000000 1.74900000 4.74006116 -1.23650100 4.72230000
0.25510000 1.33000000 1.68800000 4.86083889 -1.08044400 4.49008000
0.24900000 1.33000000 1.63100000 4.97991968 -0.89126100 4.33846000
0.24260000 1.32000000 1.57700000 5.11129431 -0.74452900 4.16328000
0.23710000 1.32000000 1.53600000 5.22986082 -0.61689600 4.05504000
0.23130000 1.30000000 1.49700000 5.36100303 -0.55100900 3.89220000
0.22620000 1.31000000 1.46000000 5.48187445 -0.41550000 3.82520000
0.22140000 1.30000000 1.42700000 5.60072267 -0.34632900 3.71020000
0.21640000 1.30000000 1.38700000 5.73012939 -0.23376900 3.60620000
0.21190000 1.30000000 1.35000000 5.85181689 -0.13250000 3.51000000
0.20730000 1.30000000 1.30400000 5.98166908 -0.01041600 3.39040000
0.20330000 1.33000000 1.27700000 6.09936055 0.13817100 3.39682000
0.19930000 1.33000000 1.25100000 6.22177622 0.20389900 3.32766000
0.19530000 1.34000000 1.22600000 6.34920635 0.29252400 3.28568000
0.19160000 1.32000000 1.20300000 6.47181628 0.29519100 3.17592000
0.18790000 1.28000000 1.18800000 6.59925492 0.22705600 3.04128000

EOF

$ECHO "done"

ELAPSED_TIME4=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME4 seconds"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-Begin Plotting with GNUPlot (Assuming That It is Installed)-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
if [ "$GP_COMMAND" = "" ]; then
    break
else
cat > gnuplot.tmp <<EOF
#!$GP_COMMAND

set term png size 1920,1080
set xrange [0:6.5]

set out 'epsi_Au.png'
set timestamp
set title 'epsi_Au'
set ylabel 'epsilon'
set xlabel 'Energy (eV)'
set yrange [-50:50]
plot "epsi_Au.dat" using 1:2 title 'epsi_x' ps 0.5, "epsi_Au.dat" using 1:3 title 'epsi_y' ps 0.5, "epsi_Au.dat" using 1:4 title 'epsi_z' ps 0.5, "Au.dat" using 4:6 ps 0.5

set out 'epsr_Au.png'
set timestamp
set title 'epsr_Au'
set ylabel 'epsilon'
set xlabel 'Energy (eV)'
set yrange [-50:50]
plot "epsr_Au.dat" using 1:2 title 'epsr_x' ps 0.5, "epsr_Au.dat" using 1:3 title 'epsr_y' ps 0.5, "epsr_Au.dat" using 1:4 title 'epsr_z' ps 0.5, "Au.dat" using 4:5 ps 0.5

set out 'eps_combined_Au.png'
set timestamp
set title 'eps_combined_Au'
set ylabel 'epsilon'
set xlabel 'Energy (eV)'
set yrange [-50:50]
plot "epsi_Au.dat" using 1:2 title 'epsi_x' ps 0.5, "epsi_Au.dat" using 1:3 title 'epsi_y' ps 0.5,  "epsi_Au.dat" using 1:4 title 'epsi_z' ps 0.5, "epsr_Au.dat" using 1:2 title 'epsr_x' ps 0.5,    "epsr_Au.dat" using 1:3 title 'epsr_y' ps 0.5, "epsr_Au.dat" using 1:4 title 'epsr_z' ps 0.5, "Au.dat" using 4:6 title 'Aui' ps 0.5, "Au.dat" using 4:5 title 'Aur' ps 0.5

set out 'eps_combined_averaged_Au.png'
set timestamp
set title 'eps_combined_averaged_Au'
set ylabel 'epsilon'
set xlabel 'Energy (eV)'
set yrange [-50:50]
plot "epsi_Au.dat" u 1:(($2+$3+$4)/3) ps 0.5, "epsr_Au.dat" u 1:(($2+$3+$4)/3) ps 0.5, "Au.dat" u 4:5 ps 0.5, "Au.dat" u 4:6 ps 0.5

EOF
$ECHO
$ECHO "  Plotting Real and Imaginary Parts of Calculated Permittivies ...\c"
$GP_COMMAND < gnuplot.tmp
$ECHO " done"
rm gnuplot.tmp
fi

ELAPSED_TIME5=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME5 seconds"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-Begin Cleaning Temp Directory-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

$ECHO "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/Au.*

$ECHO
$ECHO "$EXAMPLE_DIR: done"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:End Elapsed Time Measurement-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "It has been $ELAPSED_TIME seconds"

echo "Job Completed"

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:Au_Epsilon.sh-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-Steven E. Bopp, Materials Science & Engineering-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-September 29, 2018:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-University of California, San Diego-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-End of File-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-












