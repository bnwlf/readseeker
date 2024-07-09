# Train Data
## Data Generation

The traindata for the readseeker model where generated using an jupyter notebook named [readseeker_generate_train_test_data.ipynb](readseeker_generate_train_test_data.ipynb). The RefSeq reference data where pulled at 25.03.2024.
## Data Files

The files [train.tsv.gz](train.tsv.gz) and [dev.tsv.gz](dev.tsv.gz) represents the files used for training with the DNABERT scripts. These require the training and testdata as 6-mer token sequences. For conveniance we the all data are stored in [train_seq.tsv.gz](train_seq.tsv.gz) and [dev_seq.tsv.gz](dev_seq.tsv) as continuous sequences.
