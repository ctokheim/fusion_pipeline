out_dir=$1

# annotation files
PHOS_PATH=/liulab/ctokheim/projects/data/protein_modifications/phosphosite/3_6_2019/Phosphorylation_site_dataset 
UNIPROT_XREF=data/ptm/HUMAN_9606_idmapping.dat.gz 
SNVBOX_DIR=data/snvbox

# prepare input
python scripts/format_input/convert_to_agfusion.py -i data/fusion/tcga_fusion_calls.txt -o data/fusion/tcga_fusion_agfusion_input.txt
python scripts/format_input/convert_gene_to_transcript.py -i data/fusion/tcga_fusion_agfusion_input.txt -t data/canonical_transcripts/MANE.GRCh38.v0.9.summary.txt -o data/fusion/tcga_fusion_agfusion_input_tx_v2.txt

# run fusion annotation
cd scripts/agfusion
python cli.py batch -f ../../data/fusion/tcga_fusion_agfusion_input_tx_v2.txt -db ../../agfusion.homo_sapiens.95.db -o $out_dir/agfusion -a bellerophontes

# merge results
cd ../..
python scripts/postprocess/merge_agfusion_result.py -i $out_dir/agfusion -o $out_dir/fusion_annotation/annotated_fusions.txt

# save wt sequence info
python scripts/postprocess/save_wt_seq.py -i $out_dir/fusion_annotation/annotated_fusions.txt -o $out_dir/fusion_annotation/wt_prot_sequence.txt
python scripts/postprocess/save_wt_domains.py -d agfusion.homo_sapiens.95.db -o $out_dir/fusion_annotation/wt_prot_domains.txt

#######################
# Analyze degron motifs
#######################
mkdir -p $out_dir/motif
python scripts/analyze/regex_degron_search.py -i $out_dir/fusion_annotation/wt_prot_sequence.txt -m data/motifs/motifs_Wei_lab.txt -p $PHOS_PATH -u $UNIPROT_XREF -o $out_dir/motif/wt_motif_hits.txt
# adjust uniprot IDs to use canonical uniprot ID when possible
python scripts/postprocess/switch_uniprot_canonical.py -i $out_dir/motif/wt_motif_hits.txt -c data/ptm/uniprot_canonical_isoforms.txt -o $out_dir/motif/wt_motif_hits_canonical.txt
# create SNVBox input
python scripts/format_input/make_snvbox_input.py -i $out_dir/motif/wt_motif_hits_canonical.txt -t $SNVBOX_DIR/Transcript.txt.gz -u $SNVBOX_DIR/Uniprot_Xref.txt.gz -c $SNVBOX_DIR/CodonTable.txt.gz -o $SNVBOX_DIR/motif/wt_motif_hits_canonical.snvbox_info.txt
tail -n +2 $out_dir/motif/wt_motif_hits_canonical.snvbox_info.txt | cut -f1,6,11 > $out_dir/motif/wt_motif_hits_canonical.snvbox_info.txt 
## Manual step of geting features for motifs from snvbox

# train model for degron prediction
python scripts/degron_pred/train.py -m data/train/degron_file.txt -s data/train/degron_snvbox_annot.txt -o data/train/degron_model.pickle
# score motifs
python scripts/degron_pred/score.py -i data/train/degron_model.pickle -m $out_dir/motif/wt_motif_hits_canonical.snvbox_info.txt -s $out_dir/degron_pred/wt_motif_hits_canonical.snvbox_annot.txt -o $out_dir/degron_pred/wt_motif_hits_degron_pred.txt
# analyze impact of fusions on predicted degrons
python scripts/analyze/fusion_degron_impact.py -i $out_dir/fusion_annotation/annotated_fusions.txt -m $out_dir/motif/wt_motif_hits_canonical.txt -d $out_dir/degron_pred/wt_motif_hits_degron_pred.txt -o $out_dir/degron_pred/fusion_internal_degron_impact.txt
