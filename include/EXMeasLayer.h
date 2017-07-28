#ifndef EXMEASLAYER_H
#define EXMEASLAYER_H
//*************************************************************************
//*
//*  Measurement layer class used by RTPC12 and allows Draw() function.
//*
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TCylinder.h"
#include "EXVMeasLayer.h"
#include "KalTrackDim.h"
#include "TAttDrawable.h"

class TVTrackHit;

class EXMeasLayer : public EXVMeasLayer, public TCylinder {
public:
   // Ctors and Dtor
   EXMeasLayer(TMaterial &min,
               TMaterial &mout,
               Double_t   r0,
               Double_t   lhalf,
               Bool_t     isactive = EXVMeasLayer::kActive,
         const Char_t    *name = "RTPCML");
   virtual ~EXMeasLayer();

   // Parrent's pure virtuals that must be implemented
   virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                 const TVector3   &xv) const;
   virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
   virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
   virtual void       CalcDhDa  (const TVTrackHit &ht,
                                 const TVector3   &xv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)  const;

   // Methods for MC event generation
   virtual void       ProcessHit(const TVector3    &xx,
                                       TObjArray   &hits,
                                       bool smearing=true);

   void     SetSigmaX(Double_t v)  { fSigmaX=v; }
   void     SetSigmaZ(Double_t v)  { fSigmaZ=v; }

   inline Double_t GetSigmaX() const { return fSigmaX; }
   inline Double_t GetSigmaZ() const { return fSigmaZ; }
   
   using TObject::Draw;
   //using TAttDrawable::Draw;  //this line conflicts with a previous using declaration
   virtual void Draw(Int_t color, const Char_t *opt);

private:
   Double_t fSigmaX;   // rphi resolution
   Double_t fSigmaZ;   // z  resolution

   ClassDef(EXMeasLayer,1)
};

#endif
