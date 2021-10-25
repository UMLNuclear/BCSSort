{  
  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  long x=0;
  long n=gChain->GetEntries();
  double low, high;
  bool cflag = false;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(!!cflag && fevent->Pin1E()>0) {low = high = fevent->Pin1E(); cflag = true;}
    if(low>fevent->Pin1E()) low =   fevent->Pin1E();
    if(high<fevent->Pin1E()) high = fevent->Pin1E();
  }
} 
