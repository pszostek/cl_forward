#!/usr/bin/env bash

CC=g++ make -f Makefile clean && 
CC=g++ make -f Makefile &&
./tools/test.sh
