#! /bin/sh


time ./lagpat3d.cur

time ./traj_rev

./plot_yes.py
./set_calc.py
./select_pt.py

exit
