#ifndef EXHIT_H
#define EXHIT_H

#include "TVector3.h"
#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "EXMeasLayer.h"

class EXHit : public TVTrackHit {
public:
   // Ctors and Dtor
   EXHit(Int_t m = kMdim);
   EXHit(const EXMeasLayer &ms,
               Double_t    *x,
               Double_t    *dx,
               Double_t     b,
               Int_t        m = kMdim);
   virtual ~EXHit();

   // Parent's pure virtual functions to implement
   virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;
   virtual void       DebugPrint(Option_t *opt = "")           const;

   //By Jixie: return raw position vector
   void      SetRawXv(const TVector3 &v3) { fRawXv=v3; }
   TVector3  GetRawXv() const { return fRawXv; }
private:
	
   //By Jixie: add this vector here to store the raw (x,y,z)
   TVector3 fRawXv;

   ClassDef(EXHit,1)      // Sample hit class
};

#endif
