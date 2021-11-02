# manta-strelka-vc

A Nextflow pipeline to run [manta](https://github.com/Illumina/manta) and [strelka](https://github.com/Illumina/strelka) variant callers.

The pipeline is meant to process a pair of samples per patient: normal and tumour. Therefore, generating somatic variant calling.

## Requirements

The only requirements are to have [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) and [docker](https://docs.docker.com/get-docker/) installed in the system.

## Usage

This pipeline requires the following files:

- Input `CSV` file describing the patient ID, sample name and file paths to the input `bam` and `bai` files. The following columns are required:
    * `patient_id`: id of the patient. The file should contain two rows per patient, one for the Normal sample and the other for the tumour sample.
    * `sample_id`: id of the patient sample.
    * `tumor_normal`: an uppercase "N" or "T" to specify Normal or Tumour sample.
    * `bam`: path, ftp or s3 bucket of the `bam` mapping file.
    * `bai`: path, ftp or s3 bucket of the `bai` index file.

 An example of input `CSV` file could be found at `testdata/test_input.csv`:

```CSV
patient_id,sample_id,tumor_normal,bam,bai
pb,pb_normal,N,testdata/pb_normal.bam,testdata/pb_normal.bam.bai
pb,pb_tumor,T,testdata/pb_tumor.bam,testdata/pb_tumor.bam.bai
```

- A `genome fasta` file corresponding to the assembly used while mapping the input files.
- A `genome fasta fai` index file of the `genome fasta` provided.

### Local tests

Clone the repo:

```
git clone https://github.com/lifebit-ai/manta-strelka-vc.git
cd manta-strelka-vc
```

Then, run the pipeline:

```
nextflow main.nf --config conf/test.config
```

### CloudOS tests

Use the `--config conf/test.config` parameter and open an instance with at least 2 CPUs and 3GB of RAM.
