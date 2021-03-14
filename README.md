# VCFcontam

This module implements a method for estimating the continuous
fraction of contamination of NGS sequencing reads by parsing
relevant information from a VCF file.

## Required packages:
* [pomegranate](https://pomegranate.readthedocs.io/en/latest/)
* [pysam](https://pysam.readthedocs.io/en/latest/)
* [numpy](https://numpy.org)
* [scikit-learn](https://scikit-learn.org/stable/)
* [pandas](https://pandas.pydata.org)
```
pip3 install cython
pip3 install wheel
pip3 install pandas pysam numpy pomegranate scikit-learn
```


## Input:
* A VCF file that has been bgzipped with an accompanying tabix index file.

## Output:
* A float that represents the estimated contamination fraction.

## Serialized model:
A serialized (pickled) model is used as input into predictSVR.py. This serialized model must have
been generated using the exact save version of scikit-learn as that used in predictSVR.py.


## Usage:
The model must first be trained and the model object serialized, after which a 
VCF can be used as input to generate a contamination prediction. 

The training data can be generated as a separate step using a collection of VCFs 
with known contamination fractions using the following command:

```bash
python3 VCFcontam/generate_training_data.py --inputlist <vcflistfile> --loci <bed_file_for_feature_extraction> --feature_out <output_file_name>
```
where the file given to `--inputlist` is an unheadered two-column file where the 
first column is the (full) path to a (bgzipped and tabix-indexed) VCF file with 
a known contamination level, and the second column is a float value representing 
the known contamination level present in that VCF file. After generating the training 
data, the model can be trained using the following command:

```bash
python3 VCFcontam/trainSVR.py --training_data <feature_out_from_generate_training_data> --out <serialized_model_output_name>
```

Alternatively, the training data feature extraction and model fitting can all be 
performend in a single step by supplying `trainSVR.py` with a --vcf_metadata flag,
supplying the same two-column tab-delimited file with the paths to the VCFs and 
the known contamination level of the VCFs, like so:

```bash
python3 VCFcontam/trainSVR.py --vcf_metadata <vcflistfile> --extracted_features <outputfilefortestfeatures> --out <serialized_model_output_name> --loci <bed_file_for_feature_extraction>
```
In the above command, the extracted features are written to an outfile specified by --extracted_features, 
the serialized model output is written to a file specified by --out.

To estimate contamination from a VCF after training the model, use the following command:
```bash
python3 VCFcontam/predictSVR.py --model <trained_model_object_file> --vcf <vcf_file> --index <vcf_file_tbi_index> --loci <bed_file_for_feature_extraction>
```

## Citation:
McCartney-Melstad E, Bi K, Han J, Foo CK. 2021. VCFcontam: A Machine Learning Approach to Estimate Cross-Sample Contamination from Variant Call Data. bioRxiv doi:10.1101/2021.03.12.435007

<https://www.biorxiv.org/content/10.1101/2021.03.12.435007v1>
