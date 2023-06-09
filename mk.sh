#!/bin/sh
cd bin
rm -f wsm.x
make clean
cp Makefile_r4_g Makefile
make 
cp wsm.x ../ 
cd ..
