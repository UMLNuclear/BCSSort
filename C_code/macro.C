

TH2 *ProjectXY(TH3 *h, double low, double high){
  
  TAxis *z = h->GetZaxis();
  int binlow = z->FindBin(low);
  int binhigh = z->FindBin(high);
  z->SetRange(binlow, binhigh);
  TH2 *hist = (TH2 *)h->Project3D("xy");
  
  return hist;

}
