{

  DetHit *fhit = new DetHit;
  gChain->SetBranchAddress("DetHit", &fhit);
  long x=0;
  long n=gChain->GetEntries();
  double low, high;
  int cflag = 0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(!cflag && fhit->GetNumber()==181){
      double low = fhit->GetCharge();
      double high = fhit->GetCharge();
      cflag++;
    }
    if(fhit->GetNumber()==181){
      if(low>fhit->GetCharge()) low = fhit->GetCharge();
      if(high<fhit->GetCharge()) high = fhit->GetCharge();
    }
  }

}
