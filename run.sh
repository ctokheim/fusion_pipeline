out_dir=$1
ensembl_ver=75
mkdir -p $out_dir
mkdir -p $out_dir/agfusion
mkdir -p $out_dir/fusion_annotation
mkdir -p $out_dir/degron_pred
mkdir -p $out_dir/motif
mkdir -p $out_dir/statistical_test
mkdir -p $out_dir/statistical_test/ctype
mkdir -p $out_dir/statistical_test/degron

# annotation files
PHOS_PATH=/liulab/ctokheim/projects/data/protein_modifications/phosphosite/3_6_2019/Phosphorylation_site_dataset 
ACETYL_PATH=/liulab/ctokheim/projects/data/protein_modifications/phosphosite/3_6_2019/Acetylation_site_dataset
UBIQ_PATH=/liulab/ctokheim/projects/data/protein_modifications/phosphosite/3_6_2019/Ubiquitination_site_dataset
UNIPROT_XREF=data/ptm/HUMAN_9606_idmapping.dat.gz 
SNVBOX_DIR=data/snvbox

# prepare input
python scripts/format_input/convert_to_agfusion.py -i data/fusion/tcga_fusion_calls.txt -o data/fusion/tcga_fusion_agfusion_input.txt
python scripts/format_input/convert_gene_to_transcript.py -i data/fusion/tcga_fusion_agfusion_input.txt -t data/canonical_transcripts/MANE.GRCh38.v0.9.summary.txt -c data/canonical_transcripts/custom_transcripts.txt -o data/fusion/tcga_fusion_agfusion_input_tx_v2.txt

# run fusion annotation
cd scripts/agfusion
python cli.py batch -f ../../data/fusion/tcga_fusion_agfusion_input_tx_v2.txt -db ../../agfusion.homo_sapiens.95.db -o $out_dir/agfusion -a bellerophontes
#python cli.py batch -f ../../data/fusion/nar_tumor_fusions_agfusion_input_tx.txt -db ../../agfusion.homo_sapiens.75.db -o $out_dir/agfusion -a bellerophontes

# merge results
cd ../..
python scripts/postprocess/merge_agfusion_result.py -i $out_dir/agfusion -e $ensembl_ver -o $out_dir/fusion_annotation/annotated_fusions.txt

# save wt sequence info
python scripts/postprocess/save_wt_seq.py -i $out_dir/fusion_annotation/annotated_fusions.txt -e $ensembl_ver -o $out_dir/fusion_annotation/wt_prot_sequence.txt
python scripts/postprocess/save_wt_domains.py -d agfusion.homo_sapiens.`echo $ensembl_ver`.db -e $ensembl_ver -o $out_dir/fusion_annotation/wt_prot_domains.txt
# save fusion domain info
python scripts/postprocess/merge_agfusion_domains.py -i $out_dir/agfusion -o $out_dir/fusion_annotation/fusion_prot_domains.txt

#######################
# Analyze degron motifs
#######################
python scripts/analyze/regex_degron_search.py -i $out_dir/fusion_annotation/wt_prot_sequence.txt -m data/motifs/motifs_Wei_lab.txt -p $PHOS_PATH -a $ACETYL_PATH -u $UNIPROT_XREF -o $out_dir/motif/wt_motif_hits.txt
# adjust uniprot IDs to use canonical uniprot ID when possible
python scripts/postprocess/switch_uniprot_canonical.py -i $out_dir/motif/wt_motif_hits.txt -c data/ptm/uniprot_canonical_isoforms.txt -o $out_dir/motif/wt_motif_hits_canonical.txt
# create SNVBox input
python scripts/format_input/make_snvbox_input.py -i $out_dir/motif/wt_motif_hits_canonical.txt -t $SNVBOX_DIR/Transcript.txt.gz -u $SNVBOX_DIR/Uniprot_Xref.txt.gz -c $SNVBOX_DIR/CodonTable.txt.gz -o $out_dir/motif/wt_motif_hits_canonical.snvbox_info.txt
tail -n +2 $out_dir/motif/wt_motif_hits_canonical.snvbox_info.txt | cut -f1,6,11 > $out_dir/motif/wt_motif_hits_canonical.snvbox_input.txt 
## Manual step of geting features for motifs from snvbox

# train model for degron prediction
python scripts/degron_pred/train.py -m data/train/degron_file.txt -s data/train/degron_snvbox_annot.txt -o data/train/degron_model.pickle
# score motifs
python scripts/degron_pred/score.py -i data/train/degron_model.pickle -m $out_dir/motif/wt_motif_hits_canonical.snvbox_info.txt -s $out_dir/degron_pred/wt_motif_hits_canonical.snvbox_annot.txt -t data/train/degron_file.txt -o $out_dir/degron_pred/wt_motif_hits_degron_pred.txt
# analyze impact of fusions on predicted degrons
python scripts/analyze/fusion_degron_impact.py -i $out_dir/fusion_annotation/annotated_fusions.txt -m $out_dir/motif/wt_motif_hits_canonical.txt -d $out_dir/degron_pred/wt_motif_hits_degron_pred.txt -o $out_dir/degron_pred/fusion_internal_degron_impact.txt

######################
# Add annotations to fusion genes
######################
# -fc data/fusion/nar_tumor_fusion_calls_formatted.txt \
python scripts/postprocess/add_fusion_annotations.py \
    -i $out_dir/fusion_annotation/annotated_fusions.txt \
    -d $out_dir/degron_pred/fusion_internal_degron_impact.txt \
    -c $out_dir/degron_pred/cterm_degron_results.txt \
    -fc data/fusion/tcga_fusion_calls.txt \
    -wd $out_dir/fusion_annotation/wt_prot_domains.txt \
    -fd $out_dir/fusion_annotation/fusion_prot_domains.txt \
    -m data/mc3.v0.2.8.PUBLIC.code.filtered.small.maf \
    -cgc data/drivers/Census_allThu\ Apr\ \ 9\ 21_19_01\ 2020.tsv \
    -ot data/drivers/OG_TSG_annotation.txt   \
    -ok data/drivers/oncokb_4_3_2017.txt \
    -drug data/drugability/interactions.tsv \
    -driver data/drivers/driver_flags/PANCAN.txt \
    -u $UBIQ_PATH \
    -xref $UNIPROT_XREF \
    -o $out_dir/fusion_annotation/full_annotated_fusions.txt
#######################
# cancer type specific test
#######################
python scripts/analyze/ctype_specificity_test.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt -n 10000 --three-prime -o $out_dir/statistical_test/ctype/three_prime
python scripts/analyze/ctype_specificity_test.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt -n 10000 --three-prime --domain -o $out_dir/statistical_test/ctype/three_prime_domain
python scripts/analyze/ctype_specificity_test.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt -n 10000 -o $out_dir/statistical_test/ctype/five_prime
python scripts/analyze/ctype_specificity_test.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt -n 10000 --domain -o $out_dir/statistical_test/ctype/five_prime_domain
#######################
# Internal Degron loss
#######################
python scripts/analyze/degron_loss_test_v2.py -n 10000 -i $out_dir/fusion_annotation/full_annotated_fusions.txt --three-prime -o $out_dir/statistical_test/degron/three_prime_loss.txt
python scripts/analyze/degron_loss_test_v2.py -n 10000 -i $out_dir/fusion_annotation/full_annotated_fusions.txt --three-prime --domain -o $out_dir/statistical_test/degron/three_prime_loss_domain.txt
python scripts/analyze/degron_loss_test_v2.py -n 10000 -i $out_dir/fusion_annotation/full_annotated_fusions.txt -o $out_dir/statistical_test/degron/five_prime_loss.txt
python scripts/analyze/degron_loss_test_v2.py -n 10000 -i $out_dir/fusion_annotation/full_annotated_fusions.txt --domain -o $out_dir/statistical_test/degron/five_prime_loss_domain.txt
##############################
# Internal Degron loss partner
##############################
python scripts/analyze/degron_loss_test_partner_v2.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt --three-prime -o $out_dir/statistical_test/degron/three_prime_loss_partner.txt
python scripts/analyze/degron_loss_test_partner_v2.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt --three-prime --domain -o $out_dir/statistical_test/degron/three_prime_loss_partner_domain.txt
python scripts/analyze/degron_loss_test_partner_v2.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt -o $out_dir/statistical_test/degron/five_prime_loss_partner.txt
python scripts/analyze/degron_loss_test_partner_v2.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt --domain -o $out_dir/statistical_test/degron/five_prime_loss_partner_domain.txt
#######################
# C-terminal degron loss
#######################
python scripts/analyze/cterm_degron_loss_test.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt -m data/motifs/cterm_gps_motifs.txt -w $out_dir/fusion_annotation/wt_prot_sequence.txt -o $out_dir/statistical_test/degron/five_prime_cterm_loss.txt
python scripts/analyze/cterm_degron_loss_test.py -i $out_dir/fusion_annotation/full_annotated_fusions.txt --domain -m data/motifs/cterm_gps_motifs.txt -w $out_dir/fusion_annotation/wt_prot_sequence.txt -o $out_dir/statistical_test/degron/five_prime_cterm_loss_domain.txt
