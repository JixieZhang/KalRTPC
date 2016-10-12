#!/bin/csh -b
#this macro is to run 2000 events using several KalRTPC then compare their difference

if ($#argv < 6) then
    echo "run 2000 events using several KalRTPC then compare their difference"
    echo "uages:  exe <key> <NEvent> <Pt_min> <Pt_max> <CosTh_min> <cosTh_max> [cmp_only]"
    exit 999
endif
set cmp_only = (0)
if ($#argv >= 7) set cmp_only = ($7)

set curdir = `pwd`
cd /media/DISK500G/work/KalmanFilter/_git_Kalman/KalRTPC/tmp

set key =  "$1_1Iter" 
if ($cmp_only == 0) then
../_EXKalRTPC_1Iter_nocurveback_phicorr 0 $2 $3 $4 $5 $6; 
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(1,1,\"$key\"\)

set key =  "$1_1Iter_last" 
if ($cmp_only == 0) then
../_EXKalRTPC_1Iter_last 0 $2 $3 $4 $5 $6; 
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(0,1,\"$key\"\)

set key =  $1_2Iter_ForwardNDiagC 
if ($cmp_only == 0) then
../_EXKalRTPC_2Iter_ForwardNDiagonalC 0 $2 $3 $4 $5 $6
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(0,1,\"$key\"\)

set key =  $1_2Iter_Forward 
if ($cmp_only == 0) then
../_EXKalRTPC_2Iter_ForwardNWholeC 0 $2 $3 $4 $5 $6
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(0,1,\"$key\"\)

set key =  $1_2Iter_BackwardNDiagC 
if ($cmp_only == 0) then
../_EXKalRTPC_2Iter_BackwardNDiagonalC 0 $2 $3 $4 $5 $6
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(0,1,\"$key\"\)

set key =  $1_2Iter_Backward
if ($cmp_only == 0) then
../_EXKalRTPC_2Iter_BackwardSmoothNWholeC 0 $2 $3 $4 $5 $6
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(0,1,\"$key\"\)

set key =  $1_2Iter_Backward_nophicorr
if ($cmp_only == 0) then
../_EXKalRTPC_2Iter_BackwardSmoothNWholeC_nophicorr 0 $2 $3 $4 $5 $6
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(0,1,\"$key\"\)

cd $curdir
