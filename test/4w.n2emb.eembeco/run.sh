#!/bin/sh

if [ -z "$TESTmanbo" ]; then
   TESTmanbo="../../manbo"
fi

if [ -z "$DO_PARALLEL" ]; then
   export DO_PARALLEL=" "
fi

$DO_PARALLEL $TESTmanbo test.in < /dev/null || error

../dacdif test.log.save test.log
../dacdif test_forces.dat.save test_forces.dat
../dacdif test_positions.xyz.save test_positions.xyz
../dacdif test_properties.dat.save test_properties.dat
../dacdif test.rst.in.save test.rst.in
../dacdif test_velocities.dat.save test_velocities.dat

exit 0









