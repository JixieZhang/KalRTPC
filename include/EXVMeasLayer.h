#ifndef EXVMEASLAYER_H
#define EXVMEASLAYER_H
//*************************************************************************
//* ===================
//*  EXVMeasLayer Class
//* ===================
//*
//* The abstract measurement layer class which also inheriated from TAttDrawable.
//* TAttDrawable makes the measurement layer drawable.
//* 
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TCylinder.h"
#include "TVMeasLayer.h"
#include "TAttDrawable.h"
#include "KalTrackDim.h"
#include "TString.h"

class TVTrackHit;
class TNode;

class EXVMeasLayer : public TVMeasLayer, public TAttDrawable {
public:
  static Bool_t kActive;
  static Bool_t kDummy;

  // Ctors and Dtor
  EXVMeasLayer(TMaterial &min,
               TMaterial &mout,
               Bool_t     type = EXVMeasLayer::kActive,
               const Char_t    *name = "MeasL");
  virtual ~EXVMeasLayer();

  virtual void ProcessHit(const TVector3  &xx,
                          TObjArray &hits,
                          bool smearing=true) = 0;

  inline TString GetMLName() const { return fName; }
  inline TNode  *GetNodePtr() const { return fNodePtr; }

  inline void    SetNodePtr(TNode *nodep) { fNodePtr = nodep; }

private:
  TString  fName;      // layer name
  TNode   *fNodePtr;   // node pointer

  ClassDef(EXVMeasLayer,1)
};

#endif
