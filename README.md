
# IPIAbDiscov: A Python Package, developped at IPI, for NGS-Based Antibody Discovery
IPIAbDiscov is an open-source Python package designed to streamline the analysis of Next-Generation Sequencing (NGS) data from antibody display technologies, such as phage display or yeast display libraries. It provides a command-line interface for processing raw FASTQ files from antibody selection campaigns, enabling researchers to quantify sequence abundance, track enrichment across selection rounds, and identify promising lead candidates for therapeutic antibody development.

IPIAbDiscov allow to accelerate the transition from huge amount of NGS data to validated antibody leads, reducing reliance on traditional low-throughput screening while uncovering rare, high-potential clones often missed in conventional pipelines. Ideal for academic and biotech researchers in antibody discovery engineering.

# Key Features

FASTQ Processing: Quality trimming, filtering, and annotation of antibody sequences using tools like fastp and ANARCI for accurate numbering and germline assignment of Fab or scFv fragments.

Data Aggregation: Combine results from multiple samples or rounds to generate comprehensive repertoire tables.

Lead Selection: Automated ranking and selection of top-enriched sequences based on read counts and enrichment metrics.

Repeatability Checks: Identify and quantify repeated or duplicated sequences across datasets.


The package is lightweight, dependency-managed (via requirements.txt and Bioconda tools), and configured through YAML files and sample sheets, making it suitable for MiSeq or similar NGS runs in antibody discovery projects.
Potential Future Enhancements

Developability Assessment: Built-in filters for liabilities (e.g., glycosylation sites, cysteine residues) and integration with external databases for off-target prediction.

To further empower antibody discovery workflows, future versions could include:

Advanced Visualization: Interactive plots for repertoire diversity (e.g., Shannon entropy, clonal frequency distributions), enrichment heatmaps across selection rounds, sequence logos for CDR regions, and pairwise similarity networks.

Fold Change and Enrichment Analysis: Statistical computation of log-fold changes between pre- and post-selection rounds, with significance testing to highlight antigen-specific binders.

Machine Learning Integration: Predictive models for affinity or developability scoring, clustering of related sequences (e.g., lineage grouping), or epitope binning using sequence features.


Repository: https://github.com/bmhoan/IPIAbDiscov




# Package installtion

download abodydisco from ipi githup

#install requirements
pip install -r requirements.txt

#fastp instalation
conda install -c bioconda fastp

#anarci installation
https://github.com/oxpig/ANARCI.git
conda install -c conda-forge biopython -y
conda install -c bioconda hmmer=3.3.2 -y
cd ANARCI
python setup.py install

# How To Use AbodyDisco

cd abodydisco

python __main__.py process \
  --config config.yaml \
  --sample-sheet Miseq104/Fastq/Miseq104_SampleSheet.xlsx \
  --fastq-folder Miseq104/Fastq

python  __main__.py combine --folder Miseq104/results

python  __main__.py  pick-leads --folder Miseq104/results

python  __main__.py check-repeats --folder Miseq104/results

## Library Selection

Use one `config.yaml` for Fab, VHH, or mixed runs. The CLI accepts aliases:

```bash
python __main__.py run-all --config config.yaml --lib fab ...
python __main__.py run-all --config config.yaml --lib vhh ...
python __main__.py run-all --config config.yaml --lib mixed ...
```

New sample sheets require three minimum columns:

- `TubeBarcode`: FASTQ tube/index identifier; FASTQ names should begin with this value.
- `Sample_Name`: biological/output sample name, replacing the old `Description`
  column. Format: `target__block__round__arm__condition`.
- `library`: friendly library type, usually `fab` or `vhh`. These resolve to
  the concrete config libraries such as `standard_fab` or `vhh_full`.

Each `target + block` group should use one library. A sequencing run can mix
Fab and VHH samples, but a group such as `TargetA / Block1` cannot contain both
libraries; the validator will stop on that because clone grouping and lead
selection are library-specific.

Legacy sheets with `Sample_ID` are still accepted and converted to `TubeBarcode`
with a warning. Legacy sheets with `Description` are still accepted and
converted to `Sample_Name` with a warning. When the same target contains both
Fab and VHH samples, clone outputs are split by library, for example
`Target_Block1_standard_fab_clones.csv` and `Target_Block1_vhh_full_clones.csv`.

ML prediction is also library-aware. Fab inputs require full `HSEQ` and `LSEQ`
and can use VH+VL language-model embeddings. VHH inputs have no real light
chain: the pipeline keeps `LSEQ` empty, requires only full `HSEQ` plus `CDR3`,
and uses Delphi sequence-only models such as `--lm biophysical` or `--lm kmer`
with RF/XGBoost models when those checkpoints are registered. Successful
prediction score/label columns are merged back into each `*_final_leads.xlsx`.

# My pipeline configuration
#config.yaml - Single master configuration file

current_library: "standard_fab"   # ← CHANGE THIS LINE

general:
  base_dir: "/Users/Hoan.Nguyen/ComBio/AbodyDiscov"
  
  fastp_path: "/opt/anaconda3/bin/fastp"
  
  mafft_path: "/opt/anaconda3/bin/mafft"
  
  previous_antibodies_db: "/Users/Hoan.Nguyen/ComBio/AbodyDiscov/data/All_mAb_20251106_FACS_BLI.xlsx"
  
  output_folder: "results"



libraries:

  standard_fab:
  
    vh_barcodes:
    
      VHH-IPI-101: "TGTGCCATTTCGGGTTCCGGTGGTAGCACCTAC"
      
      H1-69: "GGTGGTATTATTCCAATTTTTGGTACTGCTAAT"
      
      H3-23_A: "TCTTATATATCTTCATCTGGTTCTACTATTTAT"
      
      H5-51: "GGTATAATCTACCCCGGTGATTCTGATACTAGA"
      
      H3-7_A: "GCTAATATCAAACAAGAAGGTTCTGAAAAGTAT"
      
      H4-34: "GGTGAAATCAATCACTCTGGTTCCACCAAC"
      
      Scaffold16: "CGTCAGGCACCGGGTAAAGAACGTGAATTAGTCAGCGCC"
      
      VOB: "AGACAAGCTCCAGGTAAGGGTCGTGAATTAGTTGCCGGT"
      
      VH_BLAM: "TACAGACAAGCTCCAGGTAAGGGTCGTGAATTAGTTGCCGGT"
      
    vl_barcodes:
    
      K1-39: "GTCACGGTATCC"
      
      K3-15: "GTGACAGTCTCG"
      
      K4-1: "GTAACCGTGTCA"
      
      K3-20: "GTTACTGTTTCTTCTGCATCTACT"
      
    vh_barcode_region: [0, -1]
    
    vl_barcode_region: [-70, -1]


  vhh_full:
  
    vh_barcodes:
    
      VHH-IPI-101: "TGTGCCATTTCGGGTTCCGGTGGTAGCACCTAC"
      
      H1-69: "GGTGGTATTATTCCAATTTTTGGTACTGCTAAT"
      
      H1-46: "GGTATTATCAAACCATCTGGTGGTTCTACTTCT"
      
      H3-23_A: "TCTTATATATCTTCATCTGGTTCTACTATTTAT"
      
      H5-51: "GGTATAATCTACCCCGGTGATTCTGATACTAGA"
      
      H5-51_A: "GGTATAATCTACCCCGGTTATTCTGATACTAGA"
      
      H3-15: "GGTCGTATTAAGAGTAAAACCGATGGTGGTACTACTGAT"
      
      H3-7_A: "GCTAATATCAAACAAGAAGGTTCTGAAAAGTAT"
      
      H4-34: "GGTGAAATCAATCACTCTGGTTCCACCAAC"
      
      H4-39: "GGTTCTATATACTATTCTGGTTCAACTTAT"
      
      DEEPASH: "ACTCGCGGCTCAACCCGCAATGGCC"
      
      DEEPASH-R: "AACAACTTTCAACAGTTTCGGCACC"
      
      Scaffold16: "CGTCAGGCACCGGGTAAAGAACGTGAATTAGTCAGCGCC"
      
      VOB: "AGACAAGCTCCAGGTAAGGGTCGTGAATTAGTTGCCGGT"
      
      VH_BLAM: "TACAGACAAGCTCCAGGTAAGGGTCGTGAATTAGTTGCCGGT"
      
    vl_barcodes:
    
      K1-39: "GTCACGGTATCC"
      
      K3-15: "GTGACAGTCTCG"
      
      K4-1_C: "GTAACCGTGTCA"
      
      K3-20: "GTTACTGTTTCTTCTGCATCTACT"
      
    vh_barcode_region: [0, -1]
    
    vl_barcode_region: [0, -1]

processing:

  yyc_nt:
  
    - "TACTACTGC"
    
    - "TATTACTGC"
    
  wgq_nt:
  
    - "TGGGGACAA"

  hcdr3_min_len: 1
  
  hcdr3_max_len: 30
  
  read_count_min: 1
  
  read_freq_min: 0.0000001   

  filters_include: []       
  
  filters_exclude: []        
  skip_processed: false

combine:   

  #Columns used for grouping in pivot table
  
  pivot_cols:
  
    - "cdr3_aa"
    
    - "vh_scaffold"
    
    - "vl_scaffold"
    
    - "cdr3_functional"
    

  #Name of the CDR3 column (used for liabilities, clustering, length filtering)
  
  cdr3_col: "cdr3_aa"

  #Filtering thresholds
  
  min_cdr3_len: 6
  
  max_cdr3_len: 30
  
  min_freq: 0.0001
  
  min_count: 5
  
  min_freq_sum: 0.0001
  
  min_reads_per_round: 5000
  
  remove_non_functional: true

  #Critical liabilities for rank adjustment
  


pick_leads:

  min_freq_first: 0.005
  
  min_freq_sum: 0.0005
  
  min_freq_table: 0.0001
  
  n_leads: 10000000
  
  critical_filtering: false
  
  priority_concentrations:
  
    - "4nM"
    
    - "20nM"
    
    - "100nM"
    
  dont_order_antigens: []

repeat_check:

  blosum_threshold: 0.8
  
  max_cdr3_len_blosum: 25
  
  n_workers: 8

fastp:

  qualified_quality_phred: 25
  
  unqualified_percent_limit: 20
  
  length_required: 50
  
  n_base_limit: 5
  
  correction: true
  
  overlap_len_require: 50
  
  overlap_diff_limit: 5
  
  thread: 8
  
  disable_adapter_trimming: true
  
  disable_trim_poly_g: true

  
#contact {Hoan.Nguyen, Andre.Teixeira}@proteininnovation.org}
