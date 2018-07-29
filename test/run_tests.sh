#!/bin/sh

# The TESTmanbo path needs to be visible from within the test folder
if [ -z "$TESTmanbo" ]; then
   TESTmanbo="../../bin/manbo"
fi

echo "*************"
echo "  Test Info"
echo "*************"
echo "  Running test with the executable:" ${TESTmanbo##*/}

if [ -z "$DO_PARALLEL" ]; then
   export DO_PARALLEL=" "
   echo "  The DO_PARALLEL variable is not set."
else
   echo "  The DO_PARALLEL variable is set as: " $DO_PARALLEL
fi

if [ -z "$OMP_NUM_THREADS" ]; then
   echo "  The OMP_NUM_THREADS variable is not set."
else
   echo "  The OMP_NUM_THREADS variable is set as: " $OMP_NUM_THREADS
fi

echo ""

# Listing tests
tests_list="4w.n2emb.eembeco.g1_bse 4w.n2emb.eembe.g1_long 4w.n3emb.eembe.g1
20w.n2emb.eembe.g1_pbc 20w.n2emb.eembe.g2_pbc 20w.n2emb.eembe.g3_pbc"

# Executing tests
for test in $tests_list
do
   echo "********************************************"
   echo "  Executing test: " $test
   echo "********************************************"
   echo "=============================================================="
   cd $test
   ./run.sh
   cd ..
done

exit 0
