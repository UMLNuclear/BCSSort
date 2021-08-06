
#include<cstdio>
#include<GCanvas.h>


GCanvas::GCanvas (Bool_t build) : TCanvas(build){
}

GCanvas::GCanvas (const char *name, const char *title, Int_t ww, Int_t wh) : TCanvas(name, title, ww, wh){}

GCanvas::GCanvas (const char *name, const char *title, Int_t form) : TCanvas(name, title, form){}

GCanvas::~GCanvas(){}

bool GCanvas::HandleKeyboard(EKeySym key) {
   bool used = false;    
   switch(key) {
      case kKey_e:
        printf("the e key was pressed!\n");
        break;
      case kKey_l:
        SetLogy(!GetLogy());
        used = true;
        break;
   }
  return used;
}

void GCanvas::HandleInput(EEventType button, Int_t x, Int_t y){
  printf("EEventType = %i; \t x = %i; \t y = %i;\n ", button,x,y);
  bool used = false;
  if(button==kKeyPress) {
    used = HandleKeyboard((EKeySym)x);
  } 
  
  if(!used) 
    TCanvas::HandleInput(button, x, y);
  else {
    if(!gPad) return;
    gPad->Modified();
    gPad->Update();
  }
}

