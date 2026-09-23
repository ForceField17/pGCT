#cat ../test.txt | awk '($52=="A[C-G]G" || $52=="A[C-T]G" || $52=="C[T-G]T" || $52=="C[T-C]T")' | cut -f 4,5,6,7,10,14,15,40 > SBSyst_distribution.txt
cat ../step5.GCT.txt | awk '($52=="A[C-G]G" || $52=="A[C-T]G" )' | cut -f 4,5,6,7,10,14,15,40 > SBSyst_distribution.txt
