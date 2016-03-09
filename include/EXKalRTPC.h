#include "TROOT.h"
#include "TApplication.h"
#include "TRint.h"
#include "TRandom.h" 
#include "TFile.h"      
#include "TTree.h"           

#define MaxHit 200

extern TFile* gMyFile;
extern TTree* gMyTree;
extern int _index_;
extern double p_rec,pt_rec,pz_rec,th_rec,ph_rec,x_rec,y_rec,z_rec;
extern double r_rec,a_rec,b_rec;
extern int npt;
extern double step_x[MaxHit],step_y[MaxHit],step_z[MaxHit];
extern double step_px[MaxHit],step_py[MaxHit],step_pz[MaxHit];
extern double step_bx[MaxHit],step_by[MaxHit],step_bz[MaxHit];
extern int step_status[MaxHit];
extern double step_x_exp[MaxHit],step_y_exp[MaxHit],step_z_exp[MaxHit];
extern double step_x_fil[MaxHit],step_y_fil[MaxHit],step_z_fil[MaxHit];

//varible from global fitter, come from g4 root tree 
extern double p0,pt0,pz0,th0,ph0,_x0_,_y0_,_z0_;
extern double p_hel,pt_hel,pz_hel,th_hel,ph_hel,x_hel,y_hel,z_hel;
extern double r_hel,a_hel,b_hel;

extern int ndf;
extern double chi2,cl;

