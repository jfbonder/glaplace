#!/bin/bash
# Final production run (one code version). Two sequential queues, one per core.
cd "$(dirname "$0")"
L=../results/v4
qA() {
  python exp01_lavrentiev_sweep.py --jobs 1 --maxn 256 > $L/exp01.log 2>&1
  python exp01_lavrentiev_sweep.py --jobs 1 --maxn 256 --big > $L/exp01_big.log 2>&1
}
qB() {
  for e in exp02_r_study exp03_bos_quadrature exp04_exponential exp05_disk_radial exp06_rough_load exp08_old_experiments exp07_double_phase_elm; do
    python $e.py > $L/$e.log 2>&1
  done
}
qA & qB & wait
