# prepare input
python scripts/format_input/convert_to_agfusion.py -i data/fusion/tcga_fusion_calls.txt -o data/fusion/tcga_fusion_agfusion_input.txt
python scripts/format_input/convert_gene_to_transcript.py -i data/fusion/tcga_fusion_agfusion_input.txt -t data/canonical_transcripts/MANE.GRCh38.v0.9.summary.txt -o data/fusion/tcga_fusion_agfusion_input_tx_v2.txt

# run fusion annotation
cd scripts/agfusion
python cli.py batch -f ../../data/fusion/tcga_fusion_agfusion_input_tx_v2.txt -db ../../agfusion.homo_sapiens.95.db -o ../../output/agfusion_tx_v2 -a bellerophontes

# merge results
cd ../..
python scripts/postprocess/merge_agfusion_result.py -i output/agfusion_tx_v2 -o output/annotated_fusions_tx_v2.txt

# save wt sequence info
python scripts/postprocess/save_wt_seq.py -i output/annotated_fusions_tx_v2.txt -o output/wt_prot_sequence.txt
python scripts/postprocess/save_wt_domains.py -d agfusion.homo_sapiens.95.db -o output/wt_prot_domains.txt
