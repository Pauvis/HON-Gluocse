# HON-Gluocse
Analysis codes for Viskaitis et al 2024 "Temporal features of blood glucose tracked by orexin neurons of behaving mice"


The workflow for 2-photon single-cell analysis :
- Cell Ca fluorescence extraction from hand drawn ROIs: "cell_f_extraction_2P_doublecheckplanes_24012022a.m"
- Align metadata with correct experiments: "ExportPreprocesseddata_2PGlusensstim08022022.m"
- Combine cells from same and nearby planes if spacialy overlapping and highly correlated: "Combine_2Pdata0802.m"
- Basic analysis and plots: "Analysis_combined_data_2P.m"
- 3D localisation of the cells: "Locate2Pcells.m"
- Align and correlate with running: "Chech_running_per_gluClasses.m"
