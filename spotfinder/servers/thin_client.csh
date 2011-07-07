#/bin/csh -f
#source /net/cci/sauter/filer_xds/setpaths.csh
source /net/cci/sauter/filer_worm/setpaths.csh
rm -rf all_results.out

foreach x (`python -c "for x in xrange(1,101): print '%03d'%x,"`)
  set file="/net/sunbird/raid1/sauter/rawdata/elspeth/p4h/P4H-semet_p4h_apo2_1_$x.img"

  # sleep length must be optimized based on server speed & number of processors
  /bin/sleep 1.60

  echo $file
  (libtbx.python -c "from spotfinder.servers.thin_client import do_main; do_main( '${file}', 'localhost', 8125 )" >> all_results.out) &

end
