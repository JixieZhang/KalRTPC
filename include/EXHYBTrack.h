#ifndef EXHYBTRACK_H
#define EXHYBTRACK_H
//*************************************************************************
//* =================
//*  EXHYBTrack Class
//* =================
//*
//* (Description)
//*   Hybrid track class for Kalman filter. Originally from example program
//*
//*************************************************************************
//
#include "TKalTrack.h"         // from KalTrackLib
#include "TAttDrawable.h"      // from Utils

class EXHYBTrack : public TKalTrack, public TAttDrawable {
public:
  EXHYBTrack(Int_t n = 1) : TKalTrack(n) {}
  ~EXHYBTrack() {}

  //using TAttDrawable::Draw;
  using TCollection::Draw;
  virtual void Draw(Int_t color, const Char_t *opt="");

  ClassDef(EXHYBTrack,1)  // Hybrid track class for Kalman Filter
};

#endif
