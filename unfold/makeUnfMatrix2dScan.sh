for etaBin in $(seq 0 2); do
    #kSample == 1 -> kPbPb
    kSample=1
    optX=77
    optY=772
    doReweight=1
    nSys=-1
    ./analysis $kSample $optX $optY $etaBin $doReweight $nSys
done


