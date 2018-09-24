#!/bin/sh

if [ -z "$TESTmanbo" ]; then
   TESTmanbo="../../bin/manbo"
fi

# make sure to set OMP_NUM_THREADS properly
if [ -z "$DO_PARALLEL" ]; then
   export DO_PARALLEL=" "
fi

$DO_PARALLEL $TESTmanbo test.in < /dev/null || error

../dacdif -r 1.0e-07 test.log.save test.log
../dacdif -r 1.0e-07 test_forces.dat.save test_forces.dat
../dacdif -r 1.0e-07 test_positions.xyz.save test_positions.xyz
../dacdif -r 1.0e-07 test_properties.dat.save test_properties.dat
../dacdif -r 1.0e-07 test.rst.in.save test.rst.in
../dacdif -r 1.0e-07 test_velocities.dat.save test_velocities.dat

rm inp*.chk

exit 0
