# peak locations for 3a
# Used a 1000kb windown around peak mind points for plotting
cp ~/group/lab/akshay/HierarBinding/nil_19022017_reAnalysis/NIL_all.peaks 3a_peaks_locations.tab

#3b
cp ~/group/lab/akshay/HierarBinding/nil_10312016_v2/nil_sequnwincer_out/Discrim_motifs_simple_reorderd.scores 3b_SeqUnwinder.tab
cp ~/group/lab/akshay/HierarBinding/nil_10312016_v2/nil_sequnwincer_out/Discrim_motifs.transfac 3b_motifs.transfac
java org.seqcode.projects.sequnwinder.utils.HeatMapMaker --vals 3b_SeqUnwinder.tab --minVal -0.4 --maxVal 0.4 --motifs 3b_SeqUnwinder_motifs.transfac --raster --nolabs

cp /storage/home/auk262/group/lab/akshay/HierarBinding/nil_19022017_reAnalysis/sequninder_MCC/MCC_Out_v2/Discrim_motifs.scores 3b_MCC.tab
cp /storage/home/auk262/group/lab/akshay/HierarBinding/nil_19022017_reAnalysis/sequninder_MCC/MCC_Out_v2/Discrim_motifs.transfac 3b_MCC_motifs.transfac
 java org.seqcode.projects.sequnwinder.utils.HeatMapMaker --vals 3b_MCC.tab --minVal -0.4 --maxVal 0.4 --motifs 3b_MCC_motifs.transfac --raster --nolabs

