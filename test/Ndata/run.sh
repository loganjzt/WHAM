mv a.out ao.out
module load gromacs-gcc
g++ -DCPLUSPLUS -Wl,-rpath -Wl,/home/ertexi/library/liblbfg/lib/ -I /home/ertexi/library/liblbfg/include/ -L /home/ertexi/library/liblbfg/lib/ -llbfgs ../uwham/load_ebWk.cpp ../uwham/save_whaminfo.cpp ../uwham/call_bfgs.cpp ../uwham/misc_fn.cpp ../uwham/getPvN.cpp ../uwham/check_PvN.cpp main.cpp
./a.out > output_Ndata.dat
