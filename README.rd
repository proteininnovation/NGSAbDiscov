

1. package installtion

download abodydisco from ipi githup

# install requirements
pip install -r requirements.txt

# fastp instalation
conda install -c bioconda fastp

#anarci installation
https://github.com/oxpig/ANARCI.git
conda install -c conda-forge biopython -y
conda install -c bioconda hmmer=3.3.2 -y
cd ANARCI
python setup.py install

2.  How To Use AbodyDisco

cd abodydisco


# VHH
nohup python /alphafold/combio/software/IPIAbDiscov/__main__.py run-all --fastq-folder /alphafold/combio/NGS/NGS_20260209_Miseq109_VHH/Fastq --sample-sheet  /alphafold/combio/NGS/NGS_20260209_Miseq109_VHH/Miseq109_SampleSheet.xlsx --config  /alphafold/combio/software/IPIAbDiscov/config_VHH.yaml  > Miseq109_VHH.log&

# FAB
nohup python /alphafold/combio/software/IPIAbDiscov/__main__.py run-all --fastq-folder /alphafold/combio/NGS/NGS_20260209_Miseq109_VHH/Fastq --sample-sheet  /alphafold/combio/NGS/NGS_20260209_Miseq109_VHH/Miseq109_SampleSheet.xlsx --config  /alphafold/combio/software/IPIAbDiscov/config.yaml  > Miseq109.log&

python __main__.py process \
  --config config.yaml \
  --sample-sheet /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/Fastq/Miseq104_SampleSheet.xlsx \
  --fastq-folder /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/Fastq

python  __main__.py combine --folder /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/results
python  __main__.py  pick-leads --folder /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/results
python  __main__.py check-repeats --folder /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/results
