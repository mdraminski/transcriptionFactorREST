#download RData files
wget -c http://zbo.ipipan.waw.pl/tools/papers/rest_2022/RData/Glioma_TCGA_Merged.RData;
wget -c http://zbo.ipipan.waw.pl/tools/papers/rest_2022/RData/TCGA_LGG_GBM_DESeq2_own_normalization_IDH_MUT_vs_WT_DESeq2_analysis_DEC_2021_plus_G4_WT_fr3.RData;
wget -c http://zbo.ipipan.waw.pl/tools/papers/rest_2022/RData/data.other.tar.gz;

tar -xzvf data.other.tar.gz;
