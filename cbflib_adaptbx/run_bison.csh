#! /bin/csh -fe
set verbose
bison cbf.stx.y -o cbf.stx.tab.c -d
mv cbf.stx.tab.c cbf_stx.c
mv cbf.stx.tab.h ../include/cbf_stx.h
