#include "EXVKalDetector.h"
#include "EXVMeasLayer.h"
#include "TVKalDetector.h"
#include "TTUBE.h"
#include "TNode.h"
#include "TRotMatrix.h"
#include "TVirtualPad.h"

//By jixie: the system think the particle carry negative charge by default
//It does not allow me to give it negative field
//Need to spend time to figure out why and solve it!

Double_t EXVKalDetector::fgBfield  = 50.;
TNode   *EXVKalDetector::fgNodePtr = 0;

ClassImp(EXVKalDetector)

EXVKalDetector::EXVKalDetector(Int_t m)
: TVKalDetector(m), fIsPowerOn(kTRUE)
{
}

EXVKalDetector::~EXVKalDetector()
{
}

TNode *EXVKalDetector::GetNodePtr()
{
  if (!fgNodePtr) {
    new TRotMatrix("rotm","rotm", 10.,80.,10.,80.,10.,80.);
    new TTUBE("Det","Det","void",6.999,7.,20.);
    fgNodePtr = new TNode("World","World","Det",0.,0.,0.,"rotm");
  }
  return fgNodePtr;
}

void EXVKalDetector::Draw(Int_t color, const Char_t *opt)
{
  if (!gPad) return;
  TNode *nodep = GetNodePtr();
  nodep->cd();
  TIter next(this);
  TObject *objp;
  while ((objp = next())) {
    TAttDrawable *dp = dynamic_cast<TAttDrawable *>(objp);
    if (dp) dp->Draw(color, opt); 
  }
  nodep->Draw("pad same");  //why use pad?
  //nodep->Draw("same");
}
