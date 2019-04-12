int ptBinning() {
  int num_reco_bins = 48;
  double dx = std::pow(10.,0.05);
  std::vector<double> reco_bins(num_reco_bins+1,0);
  reco_bins[0]=10;
  for (int i = 1; i <= num_reco_bins; ++i) {
    reco_bins[i] = reco_bins[i-1]*dx;
    cout << "i " << i << " bins " << reco_bins[i] << endl;
  }
  return 0; 
  


}
