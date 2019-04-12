#include <vector>

int testBinnning() {
  const unsigned int num_reco_bins = 48;
  double dx = std::pow(10.,0.05); //// (100/10)^(1/20) so that after 20 bins the low edge is 100, after 40 bins its 1000 etc., this number is ~1.122...
  std::vector<double> reco_bins(num_reco_bins+1,0);
  reco_bins[0]=10;
  for (unsigned int i=1; i<=num_reco_bins; i++) {
    reco_bins[i]=reco_bins[i-1]*dx;
    cout << " reco bin " << i << " pt " << reco_bins[i] << endl;
  }
  return 0;
}
