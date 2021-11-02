# manta-strelka-vc

A Nextflow pipeline to run [manta](https://github.com/Illumina/manta) and [strelka](https://github.com/Illumina/strelka) variant callers. `manta` is an open source structural variant caller developed by Illumina while `strelka` also developed by Illumina and available open source, is more small variant caller focused.

The pipeline is meant to process a pair of samples per patient: normal and tumour, so producing somatic variant calling.

## Requirements

The only requeriments are to have [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) and [docker](https://docs.docker.com/get-docker/) installed in the system.

## Usage

This pipeline requires the following files:

- Input `CSV` file describing the patient ID, sample name and file paths to the input `bam` and `bai` files. The following columns are required:
    * `patient_id`: id of the patient. Each patient should be described by two rows of this file, one for the Normal sample and the other for the tumour sample.
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

- A genome fasta file corresponding to the assembly used while mapping the input files.
- A genome fasta fai index file of the `genome fasta` provided.

### Local tests

First, you would need to download the test data prepared from the following s3 bucket:

```
s3://eu-west-1-example-data/nihr/testdata/pb_normal.bam
s3://eu-west-1-example-data/nihr/testdata/pb_normal.bam.bai
s3://eu-west-1-example-data/nihr/testdata/pb_tumor.bam
s3://eu-west-1-example-data/nihr/testdata/pb_tumor.bam.bai
s3://eu-west-1-example-data/nihr/testdata/Homo_sapiens_assembly38.fasta
s3://eu-west-1-example-data/nihr/testdata/Homo_sapiens_assembly38.fasta.fai
```

Clone the repo:

```
git clone https://github.com/lifebit-ai/manta-strelka-vc.git
cd manta-strelka-vc
```

Place the test data in the `testdata` subfolder of the repo:

```
mv ../pb_normal* testdata
mv ../pb_tumor* testdata
mv ../Homo_sapiens* testdata
```

Then, run the pipeline (important to specify the specific CPU and MEM requirements if your machine has less than the default 30 CPUs and 60GB of RAM):

```
nextflow main.nf --input testdata/test_input.csv --genome_fasta testdata/Homo_sapiens_assembly38.fasta --genome_fasta_fai testdata/Homo_sapiens_assembly38.fasta.fai --cpus_manta 2 --memory_manta 16.GB --cpus_strelka 2 --memory_strelka 16.GB -with-docker
```

### CloudOS tests

Use only the `--config conf/test.config` parameter and open an instance with at least 30 CPUs and 60GB of RAM.
