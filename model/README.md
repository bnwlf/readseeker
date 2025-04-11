---
license: apache-2.0
pipeline_tag: token-classification
---
# Model Card for Model ID

<!-- Provide a quick summary of what the model is/does. -->

This modelcard aims to be a base template for new models. It has been generated using [this raw template](https://github.com/huggingface/huggingface_hub/blob/main/src/huggingface_hub/templates/modelcard_template.md?plain=1).

## Model Details

### Model Description

<!-- Provide a longer summary of what this model is. -->



- **Developed by:** Ben Wulf
- **Model type:** DNA Transformer Model (BERT)
- **Language(s) (NLP):** DNA
- **License:** Apache 2.0
- **Finetuned from model:** ***DNABERT6*** (Yanrong Ji et al.) [@huggingface](https://huggingface.co/zhihan1996/DNA_bert_6) , [@google drive ](https://drive.google.com/file/d/1BJjqb5Dl2lNMg2warsFQ0-Xvn1xxfFXC/view?usp=sharing)

### Model Sources [optional]

<!-- Provide the basic links for the model. -->

- **Repository:** https://github.com/bnwlf/readseeker
- **Paper :** TBD


## Uses

<!-- Address questions around how the model is intended to be used, including the foreseeable users of the model and those affected by the model. -->
The ReadSeeker-Model predicts genes within read-length-sized Genomic Regions of 300BP. And therefor discriminates between coding (CDS = Label 1) and non-coding (nonCDS = Label 0) reads.

## Bias, Risks, and Limitations

The use of sequence sizes shorter than 300BP leads expectably to reduced accuracies and therefor to higher False-Positiv and False-Negative rates.


### Recommendations

<!-- This section is meant to convey recommendations with respect to the bias, risk, and technical limitations. -->

Users (both direct and downstream) should be made aware of the risks, biases and limitations of the model. More information needed for further recommendations.

## How to Get Started with the Model

The readseeker model is 100% compatible with the DNABERT python [scripts for evaluation and prediction](https://github.com/jerryji1993/DNABERT/tree/master/examples).

[More Information Needed]

## Training Details

### Training Data

The model is trained on 3 Million 300BP long Sequences randomly sampled from RefSeq Reference Sequences from viral, bacterial and mammalian origins (human, bat, pig).
1.5 million Sequences from CDS Regions and 1.5 million sequences from nonCDS regions where used for model finetuning. For more detailed information see the model 
repository at [https://github.com/bnwlf/readseeker](https://github.com/bnwlf/readseeker)


<!-- This should link to a Dataset Card, perhaps with a short stub of information on what the training data is all about as well as documentation related to data pre-processing or additional filtering. -->

[More Information Needed]

### Training Procedure

<!-- This relates heavily to the Technical Specifications. Content here should link to that section when it is relevant to the training procedure. -->

#### Preprocessing [optional]

[More Information Needed]


#### Training Hyperparameters

- **Training regime:** [More Information Needed] <!--fp32, fp16 mixed precision, bf16 mixed precision, bf16 non-mixed precision, fp16 non-mixed precision, fp8 mixed precision -->

#### Speeds, Sizes, Times [optional]

<!-- This section provides information about throughput, start/end time, checkpoint size if relevant, etc. -->

[More Information Needed]

## Evaluation

<!-- This section describes the evaluation protocols and provides the results. -->

## Citation [optional]

<!-- If there is a paper or blog post introducing the model, the APA and Bibtex information for that should go in this section. -->

**BibTeX:**

[More Information Needed]

