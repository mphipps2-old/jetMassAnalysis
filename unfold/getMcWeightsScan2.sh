for etaBin in $(seq 0 2); do
    # kSample == 1 -> kPbPb
    icent=0
    kSample=0
    nSys=-1
    ./analysis2 $kSample $icent $etaBin $nSys
done


