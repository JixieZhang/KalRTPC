#!/bin/csh -b
#this macro is to run 2000 events using several KalRTPC then compare their difference

if ($#argv < 6) then
    echo "run several version of KalRTPC then compare their difference"
    echo "uages:  exe <key> <NEvent> <Pt_min> <Pt_max> <CosTh_min> <cosTh_max> [cmp_only]"
    exit 999
endif
set cmp_only = (0)
if ($#argv >= 7) set cmp_only = ($7)

set curdir = `pwd`
cd /media/DISK500G/work/KalmanFilter/_git_Kalman/KalRTPC/testcirclefit
mkdir -p Graph

set key =  "$1_1Iter" 
if ($cmp_only == 0) then
../_EXKalRTPC_1Iter 0 $2 $3 $4 $5 $6; 
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(1,1,\"$key\"\)
root -b -q ../macros/Dependence.C+\(\"nt_$key.root\"\)
root -b -q ../DeltaXXX.C+\(\"nt_$key.root\"\)
mkdir -p $key
rename _Pr60to250.png _$key.png Graph/*.png
mv Graph/*png *png $key

set key =  $1_1Iter_half
if ($cmp_only == 0) then
../_EXKalRTPC_1Iter_half 0 $2 $3 $4 $5 $6
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(1,1,\"$key\"\)
root -b -q ../macros/Dependence.C+\(\"nt_$key.root\"\)
root -b -q ../DeltaXXX.C+\(\"nt_$key.root\"\)
mkdir -p $key
rename _Pr60to250.png _$key.png Graph/*.png
mv Graph/*png *png $key

set key =  $1_1Iter_nocurveback
if ($cmp_only == 0) then
../_EXKalRTPC_1Iter_nocurveback 0 $2 $3 $4 $5 $6
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(0,1,\"$key\"\)
root -b -q ../macros/Dependence.C+\(\"nt_$key.root\"\)
root -b -q ../DeltaXXX.C+\(\"nt_$key.root\"\)
mkdir -p $key
rename _Pr60to250.png _$key.png Graph/*.png
mv Graph/*png *png $key


set key =  $1_1Iter_nocurveback_LM
if ($cmp_only == 0) then
../_EXKalRTPC_1Iter_nocurveback_LM 0 $2 $3 $4 $5 $6
mv h.root nt_$key.root
endif
root -b -q nt_$key.root ../macros/draw.C\(0,1,\"$key\"\)
root -b -q ../macros/Dependence.C+\(\"nt_$key.root\"\)
root -b -q ../DeltaXXX.C+\(\"nt_$key.root\"\)
mkdir -p $key
rename _Pr60to250.png _$key.png Graph/*.png
mv Graph/*png *png $key


cd $curdir
