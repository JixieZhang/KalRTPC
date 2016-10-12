#!/bin/csh -fb

#key is used to tell what is the feature for all input files 
#hname# is used to tell what is the special feature of file# among all input files 
#void DrawFromMultiFiles(const char* tgStr, const char *binStr,
#	const char *cutStr="", const char* treename="t", const char *key="",
#	const char *file0=0, int color0=1, const char *hname0="0", 
#	const char *file1=0, int color1=2, const char *hname1="1", 
#	const char *file2=0, int color2=4, const char *hname2="2", 
#	const char *file3=0, int color3=6, const char *hname3="3", 
#	const char *file4=0, int color4=8, const char *hname4="4")
#	
#3 files	
root -b -q 'DrawFromMultiFiles.C+("pt0-pt_rec","(100,-0.03,0.07)","","t",  "dPt_for_Pt50to70_All", "nt_50to70MeV_All_1Iter.root",1, "1Iter", "nt_50to70MeV_All_2Iter_Forward.root",2, "2Iter_Forward", "nt_50to70MeV_All_2Iter_Backward.root",4, "2Iter_Backword")'

root -b -q 'DrawFromMultiFiles.C+("th0-th_rec","(100,-0.2,0.2)","","t", "dTheta_for_Pt50to70_All", "nt_50to70MeV_All_1Iter.root",1, "1Iter", "nt_50to70MeV_All_2Iter_Forward.root",2, "2Iter_Forward", "nt_50to70MeV_All_2Iter_Backward.root",4, "2Iter_Backword")'

root -b -q 'DrawFromMultiFiles.C+("ph0-ph_rec","(100,-0.5,0.5)","","t",   "dPhi_for_Pt50to70_All", "nt_50to70MeV_All_1Iter.root",1, "1Iter", "nt_50to70MeV_All_2Iter_Forward.root",2, "2Iter_Forward", "nt_50to70MeV_All_2Iter_Backward.root",4, "2Iter_Backword")'



#5 files	
root -b -q 'DrawFromMultiFiles.C+("pt0-pt_rec","(100,-0.03,0.07)","","t",  "dPt_for_Pt50to70_ThAll", "nt_50to70MeV_All_1Iter.root",1, "1Iter", "nt_50to70MeV_All_2Iter_Forward.root",2, "2Iter_Forward", "nt_50to70MeV_All_2Iter_ForwardNDiagC.root",4, "2Iter_ForwardNDiagC", "nt_50to70MeV_All_2Iter_Backward.root",6, "2Iter_Backword", "nt_50to70MeV_All_2Iter_BackwardNDiagC.root",8, "2Iter_BackwordNDiagC")'

root -b -q 'DrawFromMultiFiles.C+("th0-th_rec","(100,-0.2,0.2)","","t", "dTheta_for_Pt50to70_ThAll", "nt_50to70MeV_All_1Iter.root",1, "1Iter", "nt_50to70MeV_All_2Iter_Forward.root",2, "2Iter_Forward", "nt_50to70MeV_All_2Iter_ForwardNDiagC.root",4, "2Iter_ForwardNDiagC", "nt_50to70MeV_All_2Iter_Backward.root",6, "2Iter_Backword", "nt_50to70MeV_All_2Iter_BackwardNDiagC.root",8, "2Iter_BackwordNDiagC")'

root -b -q 'DrawFromMultiFiles.C+("ph0-ph_rec","(100,-0.5,0.5)","","t",   "dPhi_for_Pt50to70_ThAll", "nt_50to70MeV_All_1Iter.root",1, "1Iter", "nt_50to70MeV_All_2Iter_Forward.root",2, "2Iter_Forward", "nt_50to70MeV_All_2Iter_ForwardNDiagC.root",4, "2Iter_ForwardNDiagC", "nt_50to70MeV_All_2Iter_Backward.root",6, "2Iter_Backword", "nt_50to70MeV_All_2Iter_BackwardNDiagC.root",8, "2Iter_BackwordNDiagC")'
