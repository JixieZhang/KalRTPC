# KalRTPC
Kalman Filter program for RTPC12
This project is built to test how well Kalman Filter can do 
in fitting RTPC tracks.
To get this source code: 
git clone git@github.com:jixie/KalRTPC.git 

This program depends on 3 other software packages:
1) root, you can download it from root.cern.ch, any version should be fine.
2) clhep, I used 2.1.0.1 compiled with vc10 in windows xp, and 2.1.3.1 in CentOS_6.5
3) KalmanFilter libraries originally from KEK, modified by me. 
You can also get it from github:  
git clone git@github.com:jixie/KalmanFilter.git 

Once you download the source code, type 'source setup.csh' to define where you want to 
install KalmanFilter package, then type 'make;make install' to install it.

To run this in windows:

vc10 project files are also aviable. You have to prepare g77,
clhep and root by youselves. Then put KalmanFilter and KalRTPC into 
the same directory. Then open the solution file with VC10. 
