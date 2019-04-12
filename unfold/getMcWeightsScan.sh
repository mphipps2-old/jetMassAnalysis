for icent in $(seq 0 6); do
#for icent in $(seq 0 0); do
#    for etaBin in $(seq 0 2); do
	# kSample == 1 -> kPbPb
#	icent=1
	kSample=1
	nSys=-1
      	./analysis $kSample $icent $nSys
#	./analysis $kSample $icent $etaBin $nSys
 #   done
done

