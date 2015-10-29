#!/usr/bin/env bash

CC=g++ make -f ps_makefile clean && 
CC=g++ make -f ps_makefile &&
./tools/test.sh
