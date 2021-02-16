TEST_SCRIPT=$PWD/cctbx_do_test.sh
TEST_OUT=$PWD/cctbx_test_result

unset CCTBX_PYTHON
CCTBX_PYTHON=$(which cctbx.python)
if [ -z $CCTBX_PYTHON ]
then
  echo
  echo "Please activate your cctbx environment and try again"
  return
fi



CCTBX_ROOT=${CCTBX_PYTHON%/build/bin/cctbx.python}

echo "source $CCTBX_ROOT/build/conda_setpaths.sh" > $TEST_SCRIPT
echo "mpirun libtbx.python -c \"import socket; print(socket.gethostname())\" " >> $TEST_SCRIPT

echo > $TEST_OUT

echo bsub -n 24 -q psdebugq -o $TEST_OUT "source $TEST_SCRIPT"
bsub -n 24 -q psdebugq -o $TEST_OUT "source $TEST_SCRIPT"

while [[ $(wc -l $TEST_OUT | awk '{print $1}') -le 1 ]]
do
  echo
  echo "Waiting for test to run"
  sleep 10
done

echo
echo "Test running. Waiting 60 seconds to finish."
sleep 60

if [[ $(grep "^psana" $TEST_OUT | wc -l) -eq 24 ]]
then
  echo "OK"
  rm $TEST_SCRIPT $TEST_OUT
else
  echo "Test failed. Output in $TEST_OUT"
fi

unset TEST_SCRIPT TEST_OUT
