#!/bin/bash

# change this if needed
#EXAT='../../src/exat.py'

if [ "$1" == "verbose" ]; then
  outstream='/dev/stdout'
else
  outstream='/dev/null'
fi

# Diffing function
#
diffing()
{
 tst=$1
 ref=$2
 # Always diff test ref (ignoring # comments) 
 diff -I '^#' -I '^Generated' $1 $2 > diff_${1}_${2}
 result=$?
 if test $result -eq 0; then
   rm diff_${1}_${2}
 fi
 echo $result
}

# Test function 
difftest()
{
  ntest=$1; shift
  nerr=0
  diffs=()
  echo ' ==================================================== '
  echo " Diffing results for test ${ntest}: "
  for arg in $@; do
    tstfile="test${ntest}.${arg}"
    reffile="ref${ntest}.${arg}"
    diffile="diff_${tstfile}_${reffile}"
    echo " ... diffing $tstfile $reffile"
    res=`diffing $tstfile $reffile`
    if test $res -ne 0; then
      nerr=`expr $nerr + 1`
      diffs="$diffs \n $diffile" 
    fi
  done
  # collect errors
  if test $nerr -gt 0; then 
    echo ' POSSIBLE FAILURE ... check these files: ' 
    echo -e $diffs
    echo
    result=1
  else
    echo ' PASSED'
    result=0
  fi
  echo ' ==================================================== '
# sleep 1
  return $result
}



##MAIN

# Set gdvh36 version
#EXAT="$EXAT -g g16"

## Run all tests for gdvh36
ntests=0
nerrors=0


echo ' ==================================================== '
echo '                      G16  TESTS                      '
echo ' ==================================================== '
echo -e  " Running test 0: "
# No seltran, no other things
$EXAT -v -i NMR4_Vac.log -o test0 | tee test0.EXAT.out > $outstream
echo
difftest 0 results.out diag.dat matrix.dat
a=$? ; nerrors=`expr $nerrors + $a`
ntests=`expr $ntests + 1`

echo " Running test 1: "
# Seltran
$EXAT -v Test1.in -i NMR4_Vac.log --seltran -o test1 | tee test1.EXAT.out > $outstream
echo
difftest 1 results.out diag.dat matrix.dat
a=$? ; nerrors=`expr $nerrors + $a`
ntests=`expr $ntests + 1`

echo " Running test 2: "
# Seltran, selecting only some chromophores 
$EXAT -v Test2.in -i NMR4_Vac.log --seltran -o test2 | tee test2.EXAT.out > $outstream
echo
difftest 2 results.out diag.dat matrix.dat
a=$? ; nerrors=`expr $nerrors + $a`
ntests=`expr $ntests + 1`

echo " Running test 3: "
# Anadipo + geometry
$EXAT -v Test3.in -i NMR4_Vac.log --seltran \
 --anadipo Test3.anadipo.in -o test3 | tee test3.EXAT.out > $outstream
echo
difftest 3 results.out diag.dat matrix.dat dipo.out geometry.xyz
a=$? ; nerrors=`expr $nerrors + $a`
ntests=`expr $ntests + 1`

echo " Running test 4: "
# Reorient + selecting chromophores
$EXAT -v Test4.in -i NMR4_Vac.log --seltran --reorient 1 2 --mag -o test4 | tee test4.EXAT.out > $outstream
echo
difftest 4 results.out diag.dat matrix.dat dipo.out magdipo.out
a=$? ; nerrors=`expr $nerrors + $a`
ntests=`expr $ntests + 1`


echo
echo ' Test completed'
echo ""
printf " %3d tests completed  \n" $ntests
printf " %3d tests failed \n"  $nerrors


echo 
echo -e 'END\n'

exit 0 

