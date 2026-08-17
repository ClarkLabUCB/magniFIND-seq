# Trusted configuration for the uncorrected PBMC/K562-r joint embedding.
INPUT_MANIFEST="embedding_inputs.tsv"
FUSION_CANDIDATES="../../02_fusion_analysis/output/fusion_detection/fusion_scan/magnifind_seq_BCR_ABL1_expected_orientation_summary_all_samples.tsv"
OUT_DIR="../output/pbmc_k562r_embedding"
K562_R_GROUP="K562-R"
SEED=20260730
PCA_COMPONENTS=30
JOINT_DIMS=20
K562_DIMS=15
JOINT_VARIABLE_FEATURES=3000
K562_VARIABLE_FEATURES=2000
