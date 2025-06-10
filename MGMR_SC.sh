#  Define some importants shortcuts for running the MGMR codes 
#@echo off is in LINUX:  scriptname > /dev/null
#
# Make sure the following line is set correctly in .bashrc  in your top folder to point to the MGMR folder on your system
# As an example:       export LL_Base=/home/olaf/MGMR3D
# 
#
# Test if the LOFLI system was installed already
if [[ ! -z "${M_ProgDir}" ]]; then
   echo MGMR system variables have already been defined
else
   echo Define MGMR system variables and PATH extensions
  
   ###### User installation dependent settings  ###############################
   export M_ProgDir="${Home}/MGMR3D/Program"
   export M_GLEscr="${Home}/MGMR3D/Program"
   export FFTlib="-lm ${Home}/NumLib/bin/libfftpack5.1d.a"
   export BLASlib="/usr/lib64/libopenblas.so"
   export AntennaFun="${Home}/LOFLI/AntenFunct/v2-"
   
   ###### User installation independent  ###############################
   export FCFLAGS="-ggdb -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -finit-real=inf"
   echo 'Program folder=' ${M_ProgDir}
   echo 'AntennaFunction=' ${AntennaFun}

fi
export RunFolder=$(pwd)
export Workdir=${RunFolder}
echo 'working directory=' ${RunFolder}

