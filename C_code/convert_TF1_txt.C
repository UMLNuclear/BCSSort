void convert_TF1_txt(TF1 *fx, double upper){

  std::ofstream ofile;
  ofile.open(Form("%s.dat",fx->GetName()));
  for(double i=0;i<=upper;i=i+0.1){
    ofile << i << "\t" << fx->Eval(i) << std::endl;
  }

  ofile.close();

}
