#!/bin/bash

rsync -avz --exclude-from=${HOME}/research/etc/rsync-excludes \
    -e ssh dk@hpc.msu.edu:src/self-assembly/expr/$1/$2 var/$1/
