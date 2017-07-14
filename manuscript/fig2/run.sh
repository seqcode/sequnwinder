# 2b(right panel)
cp ~/group/lab/akshay/HierarBinding/simulation_v7/specific_simulations/simulateOverlap_sequnwinder_04212017/Discrim_motifs_simple.scores ./
cp ~/group/lab/akshay/HierarBinding/simulation_v7/specific_simulations/simulateOverlap_sequnwinder_04212017/Discrim_motifs.transfac 2b_SeqUnwinder_motifs.transfac
java org.seqcode.projects.sequnwinder.utils.HeatMapMaker --vals 2b_SeqUnwinder.tab --minVal -0.4 --maxVal 0.4 --motifs 2b_SeqUnwinder_motifs.transfac --raster --nolabs

# 2b(left panel)
cp ~/group/lab/akshay/HierarBinding/simulation_v7/specific_simulations/simulateOverlap_mcc_04212017/Discrim_motifs_simple.scores ./
cp ~/group/lab/akshay/HierarBinding/simulation_v7/specific_simulations/simulateOverlap_mcc_04212017/Discrim_motifs.transfac 2b_MCC.transfac 
grep "A_" Discrim_motifs_simple.scores > tmp
grep "B_" Discrim_motifs_simple.scores >> tmp
grep "C_" Discrim_motifs_simple.scores >> tmp
grep "X_" Discrim_motifs_simple.scores >> tmp
grep "Y_" Discrim_motifs_simple.scores >> tmp
awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$2"\t"$6}' tmp > 2b_MCC.tab
rm tmp
rm Discrim_motifs_simple.scores
java org.seqcode.projects.sequnwinder.utils.HeatMapMaker --vals 2b_MCC.tab --minVal -0.4 --maxVal 0.4 --motifs 2b_MCC.transfac --raster --nolabs


# 2d
cp ~/group/lab/akshay/HierarBinding/simulation_v7/plots/toPlot.mat F1_scores.tab
cp ~/group/lab/akshay/HierarBinding/simulation_v7/plots/plot.R F1_toplot.R
Rscript F1_toplot.R

#2e and 2f

