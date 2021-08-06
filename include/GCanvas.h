#ifndef __GCanvas_H__
#define __GCanvas_H__

#include <TCanvas.h>

#include<KeySymbols.h>

class GCanvas : public TCanvas {

  public:
    GCanvas (Bool_t build=kTRUE);
    GCanvas (const char *name, const char *title, Int_t ww, Int_t wh);
    GCanvas (const char *name, const char *title="", Int_t form=1);  
    virtual ~GCanvas();

    virtual void HandleInput(EEventType button, Int_t x, Int_t y);
    bool HandleKeyboard(EKeySym key);

  private:

    ClassDef(GCanvas,0)
};

#endif
