import sys
import tensorflow as tf
# Limit Cores to 50 to keep it analog to FragGenScan
tf.config.threading.set_intra_op_parallelism_threads(50)
tf.config.threading.set_inter_op_parallelism_threads(50)

import keras
import pandas as pd
from pathlib import Path


from genomic_benchmarks.models.tf import vectorize_layer





model, pathtodata, ref, sample = sys.argv[-4:]
pathtodata = Path( pathtodata)

#DATASET = 'demo_coding_vs_intergenomic_seqs'
VERSION = 0
BATCH_SIZE = 64
EPOCHS = 10

CLASSES = ['intergenomic_seqs', 'coding_seqs']
NUM_CLASSES = len(CLASSES)

vectorize_layer.set_vocabulary(['', '[UNK]', 'a', 't', 'g', 'c'])
VOCAB_SIZE = len(vectorize_layer.get_vocabulary())


def vectorize_text(text, label):
  text = tf.expand_dims(text, -1)
  return vectorize_layer(text)-2, label


labels=list()
SEQ_PATH = pathtodata / "test"
model =  keras.saving.load_model(model)
test_dset = tf.keras.preprocessing.text_dataset_from_directory(
        SEQ_PATH ,
        batch_size=BATCH_SIZE,
        class_names=CLASSES,
        shuffle=False)

test_ds =  test_dset.map(vectorize_text)
for features,lab in test_ds:
    labels += lab.numpy().tolist()
# Apply tf.sigmoid to mimic keras model evaluate
raw_prediction = tf.sigmoid(model.predict(test_ds)).numpy().flatten()
prediction = (raw_prediction >= 0.5).astype("int32").tolist()
dc_pd = pd.DataFrame({"Reference":ref,"Sample":sample,"Prediction":raw_prediction.tolist(),"Label":labels,"PredictedLabel":prediction})
dc_pd.to_csv(f"simple_model_results/smr_{ref}_{sample}.tsv",sep="\t",index=False)