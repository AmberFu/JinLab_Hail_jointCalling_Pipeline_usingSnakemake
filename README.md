# JinLab_Hail_jointCalling_Pipeline_usingSnakemake
TEST Snakemake Pipeline for Hail Joint-Calling (run_combiner()) Pipeline

Follow [snakemake - Best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html)

### Submit jobs through SBUB

```
bsub -q general -G compute-jin810 \
-J hail_jointcall -N -u fup@wustl.edu \
-R 'affinity[core(5)] span[ptile=6] rusage[mem=25GB]' \
-oo master_out.txt -eo master_err.txt \
-g /fup/jobGroup_2_snakemake \
-a 'docker(spashleyfu/ubuntu20_snakemake:bamMetrics)' \
snakemake --use-conda -j 50 \
--printshellcmds --show-failed-logs \
--cluster-config $PWD/config/cluster.json \
--cluster "bsub -q general -G compute-jin810 -Ne -o {cluster.log} -e {cluster.err} -M {cluster.mem} -n {cluster.core} -R {cluster.resources} -g {cluster.jobgroup} -a 'docker({cluster.image})'" \
-s $PWD/workflow/snakefile
```

### Cluster setting

> Updated `cluster.json`: https://github.com/AmberFu/JinLab_Hail_jointCalling_Pipeline_usingSnakemake/blob/main/config/cluster.json

```json
{
    "__default__" :
    {
        "core" : 2,
        "mem" : "10GB",
        "resources" : "\"rusage[mem=10GB]\"",
        "image" : "spashleyfu/ubuntu20_snakemake:bamMetrics",
        "log" : "logs/out.txt",
        "err" : "logs/err.txt",
        "jobgroup": "/fup/jobGroup_50_snakemake"
    },
    "gvcf2gvcfbgz" : {
        "core" : 6,
        "image" : "spashleyfu/ubuntu20_snakemake:samtools",
        "log" : "logs/gvcf2gvcfbgz_out.txt",
        "err" : "logs/gvcf2gvcfbgz_err.txt"
    },
    "generate_tabix" : {
        "core" : 6,
        "image" : "spashleyfu/ubuntu20_snakemake:samtools",
        "log" : "logs/generate_tabix_out.txt",
        "err" : "logs/generate_tabix_err.txt"
    },
    "hail_run_combiner" : {
        "core" : 32,
        "mem" : "128GB",
        "resources" : "\"rusage[mem=128GB]\"",
        "image" : "spashleyfu/ubuntu18_vep104:hail_snakemake",                                                                                 
        "log" : "logs/hail_run_combiner_out.txt",
        "err" : "logs/hail_run_combiner_err.txt"
    }
}
```


### Detail:

#### Outline:

* [Submit ALL chromosome](#Submit-ALL-chromosome)
* [Test submit jobs for chr1](#test-submit-jobs-for-chr1)
* [How to BGZ a GVCF file and TABIX the file](#how-to-bgz-a-gvcf-file-and-tabix-the-file)
* [How to run Hail Joint-Calling using BASH and BSUB](#how-to-run-hail-joint-calling-using-bash-and-bsub)


#### Submit ALL chromosome:

> Date: 021/11/10
> 
> PID: 951945

```
// dry run
(base) fup@compute1-exec-135:/storage1/fs1/jin810/Active/Neuropathy_WGS_2021May/JinLab_Hail_jointCalling_Pipeline_usingSnakemake$ snakemake -n
...
Job counts:
	count	jobs
	1	all
	8419	generate_tabix
	8418	gvcf2gvcfbgz
	24	hail_run_combiner
	16862


// SUBMIT
[fup@compute1-client-3 JinLab_Hail_jointCalling_Pipeline_usingSnakemake]$ bsub -q general -G compute-jin810 \
> -J hail_jointcall -N -u fup@wustl.edu \
> -R 'affinity[core(5)] span[ptile=6] rusage[mem=25GB]' \
> -oo master_out.txt -eo master_err.txt \
> -g /fup/jobGroup_2_snakemake \
> -a 'docker(spashleyfu/ubuntu20_snakemake:bamMetrics)' \
> snakemake --use-conda -j 50 \
> --printshellcmds --show-failed-logs \
> --cluster-config $PWD/config/cluster.json \
> --cluster "bsub -q general -G compute-jin810 -Ne -o {cluster.log} -e {cluster.err} -M {cluster.mem} -n {cluster.core} -R {cluster.resources} -g {cluster.jobgroup} -a 'docker({cluster.image})'" \
> -s $PWD/workflow/snakefile
Job <951945> is submitted to queue <general>.
```

#### Test submit jobs for chr1

> Date: 2021/11/10
> 
> PID: 949643

```
[fup@compute1-client-3 JinLab_Hail_jointCalling_Pipeline_usingSnakemake]$ bsub -q general -G compute-jin810 \
> -J hail_jointcall -N -u fup@wustl.edu \
> -R 'affinity[core(5)] span[ptile=6] rusage[mem=25GB]' \
> -oo out.txt -eo err.txt \
> -g /fup/jobGroup_2_snakemake \
> -a 'docker(spashleyfu/ubuntu20_snakemake:bamMetrics)' \
> snakemake --use-conda -j 30 --show-failed-logs \
> --keep-going --printshellcmds \
> --cluster-config $PWD/config/cluster.json \
> --cluster "bsub -q general -G compute-jin810 -o {cluster.log} -e {cluster.err} -M {cluster.mem} -n {cluster.core} -R {cluster.resources} -g {cluster.jobgroup} -a 'docker({cluster.image})'" \
> -s $PWD/workflow/snakefile
Job <949643> is submitted to queue <general>.

// snakemake dry run result:
Job counts:
	count	jobs
	1	all
	366	generate_tabix
	366	gvcf2gvcfbgz
	1	hail_run_combiner
	734
```

Snakemake Rule for BGZIP VCF and create TABIX for VCF: `workflow/rules/bgzip_tabix.smk`

```python
rule gvcf2gvcfbgz:
    input: config['subset_gvcf_dir'] + "/chr{interval}/{sample}_{interval}_germline.g.vcf"
    output: config['subset_gvcf_dir'] + "/chr{interval}/{sample}_{interval}_germline.g.vcf.bgz"
    threads: config['gvcf2gvcfbgz']['threads']
    resources:
        mem_mb = config['gvcf2gvcfbgz']['mem_mb']
    shell: "bgzip --threads {threads} -c {input} > {output}"

rule generate_tabix:
    input: config['subset_gvcf_dir'] + "/chr{interval}/{sample}_{interval}_germline.g.vcf.bgz"
    output: config['subset_gvcf_dir'] + "/chr{interval}/{sample}_{interval}_germline.g.vcf.bgz.tbi"
    threads: config['generate_tabix']['threads']
    resources:
        mem_mb = config['generate_tabix']['mem_mb']
    shell: "tabix -p vcf -f {input}"
```

[[back to outline]](#outline)

Snakemake Rule for Hail Joint-Calling: `workflow/rules/hail_jointcall.smk`

```python
rule hail_run_combiner:
    input: 
        gvcf_bgz = expand(config['subset_gvcf_dir'] + "/chr{{interval}}/{sample}_{{interval}}_germline.g.vcf.bgz", sample=SAMPLES),
        gvcf_tbi = expand(config['subset_gvcf_dir'] + "/chr{{interval}}/{sample}_{{interval}}_germline.g.vcf.bgz.tbi", sample=SAMPLES),
    output: directory(config['project_outdir'] + "/Neuropathy_batch1_chr{interval}_runCombiner.mt")
    threads: config['hail_run_combiner']['threads']
    resources:
        mem_mb = config['hail_run_combiner']['mem_mb']
    params:
        spark_mem = config['hail_run_combiner']['spark_mem'],
        output_temp = config['project_outdir'] + "/hail_temp_chr{interval}"
    script: "../scripts/hail_jointCalling.py"
```

Python script for Hail Joint-Calling: `workflow/scripts/hail_jointCalling.py`

```python
import os
import time
import pyspark
import hail as hl

start_time = time.time()

threads = snakemake.threads - 2
memory = snakemake.params.spark_mem
hail_jars = "/opt/conda/lib/python3.7/site-packages/hail/backend/hail-all-spark.jar"
conf = pyspark.SparkConf().setAll([
    ('spark.master', 'local[{}]'.format(threads)),
    ('spark.app.name', 'Hail'),
    ('spark.jars', str(hail_jars)),
    ('spark.driver.extraClassPath', str(hail_jars)),
    ('spark.executor.extraClassPath', './hail-all-spark.jar'),
    ('spark.serializer', 'org.apache.spark.serializer.KryoSerializer'),
    ('spark.kryo.registrator', 'is.hail.kryo.HailKryoRegistrator'),
    ### https://discuss.hail.is/t/turning-run-combiner-performance-for-hail-local-mode/2318
    ('spark.driver.memory', str(memory)),
    ('spark.executor.memory', str(memory)),
    ])

### Using sc:
sc = pyspark.SparkContext(conf=conf)
hl.init(default_reference='GRCh38',sc=sc)

### Input: (list of files)
inputs = snakemake.input.gvcf_bgz
### Output:
temp_folder = str(snakemake.params.output_temp)
#if os.path.exists(temp_folder) and os.path.isdir(temp_folder):
#    os.system("rm -rf {}".format(temp_folder)):
os.makedirs(temp_folder, mode=777, exist_ok=True)
output_mt = snakemake.output

hl.experimental.run_combiner(inputs, 
                             use_genome_default_intervals=True,
                             out_file=output_mt, 
                             tmp_path=temp_folder,
                             overwrite=True,
                             reference_genome='GRCh38')

processing_time_in_seconds = time.time() - start_time
print("--- %s minute ---" % (processing_time_in_seconds/60))
```

[[back to outline]](#outline)


#### How to BGZ a GVCF file and TABIX the file

```
(base) fup@compute1-exec-283:/storage1/fs1/jin810/Active/Neuropathy_WGS_2021May$ time bgzip TWHJ-PNRR-10001-10001_1_germline.g.vcf && tabix -p vcf TWHJ-PNRR-10001-10001_1_germline.g.vcf.gz
real	1m40.566s
user	1m35.444s
sys	0m4.344s

// Test performance:
// BGZIP
(base) fup@compute1-exec-136:/storage1/fs1/jin810/Active/Neuropathy_WGS_2021May/JinLab_Hail_jointCalling_Pipeline_usingSnakemake$ bgzip -h

Version: 1.9
Usage:   bgzip [OPTIONS] [FILE] ...
Options:
   -b, --offset INT           decompress at virtual file pointer (0-based uncompressed offset)
   -c, --stdout               write on standard output, keep original files unchanged
   -d, --decompress           decompress
   -f, --force                overwrite files without asking
   -h, --help                 give this help
   -i, --index                compress and create BGZF index
   -I, --index-name FILE      name of BGZF index file [file.gz.gzi]
   -l, --compress-level INT   Compression level to use when compressing; 0 to 9, or -1 for default [-1]
   -r, --reindex              (re)index compressed file
   -g, --rebgzip              use an index file to bgzip a file
   -s, --size INT             decompress INT bytes (uncompressed size)
   -@, --threads INT          number of compression threads to use [1]
   -t, --test                 test integrity of compressed file

(base) fup@compute1-exec-136:/storage1/fs1/jin810/Active/Neuropathy_WGS_2021May/JinLab_Hail_jointCalling_Pipeline_usingSnakemake$ time bgzip -c TWHJ-PNRR-10001-10001_1_germline.g.vcf > TWHJ-PNRR-10001-10001_1_germline.g.vcf.bgz

real	1m45.684s
user	1m40.014s
sys	0m4.066s
(base) fup@compute1-exec-136:/storage1/fs1/jin810/Active/Neuropathy_WGS_2021May/JinLab_Hail_jointCalling_Pipeline_usingSnakemake$ time bgzip --threads 6 -c TWHJ-PNRR-10001-10001_1_germline.g.vcf > TWHJ-PNRR-10001-10001_1_germline.g.vcf.bgz

real	0m38.142s
user	1m45.398s
sys	0m5.288s

// TABIX
(base) fup@compute1-exec-136:/storage1/fs1/jin810/Active/Neuropathy_WGS_2021May/JinLab_Hail_jointCalling_Pipeline_usingSnakemake$ tabix -h

Version: 1.9
Usage:   tabix [OPTIONS] [FILE] [REGION [...]]

Indexing Options:
   -0, --zero-based           coordinates are zero-based
   -b, --begin INT            column number for region start [4]
   -c, --comment CHAR         skip comment lines starting with CHAR [null]
   -C, --csi                  generate CSI index for VCF (default is TBI)
   -e, --end INT              column number for region end (if no end, set INT to -b) [5]
   -f, --force                overwrite existing index without asking
   -m, --min-shift INT        set minimal interval size for CSI indices to 2^INT [14]
   -p, --preset STR           gff, bed, sam, vcf
   -s, --sequence INT         column number for sequence names (suppressed by -p) [1]
   -S, --skip-lines INT       skip first INT lines [0]

Querying and other options:
   -h, --print-header         print also the header lines
   -H, --only-header          print only the header lines
   -l, --list-chroms          list chromosome names
   -r, --reheader FILE        replace the header with the content of FILE
   -R, --regions FILE         restrict to regions listed in the file
   -T, --targets FILE         similar to -R but streams rather than index-jumps

(base) fup@compute1-exec-136:/storage1/fs1/jin810/Active/Neuropathy_WGS_2021May/JinLab_Hail_jointCalling_Pipeline_usingSnakemake$ time tabix -p vcf -f TWHJ-PNRR-10001-10001_1_germline.g.vcf.bgz 

real	0m34.631s
user	0m28.649s
sys	0m0.281s
```
[[back to outline]](#outline)

#### How to run Hail Joint-Calling using BASH and BSUB

BSUB file

```bash
#!/bin/bash
#BSUB -n 24
#BSUB -R 'rusage[mem=256GB]'
#BSUB -q general
#BSUB -G compute-jin810
#BSUB -J hail_02_combiner
#BSUB -N
#BSUB -u fup@wustl.edu
#BSUB -a 'docker(spashleyfu/ubuntu18_vep104:hail_gsutil)'
#BSUB -oo /storage1/fs1/jin810/Active/fup/test_pyspark_02/LOGs/test_hail_combiner_local.out
#BSUB -eo /storage1/fs1/jin810/Active/fup/test_pyspark_02/LOGs/test_hail_combiner_local.err

HAIL_HOME=$(pip3 show hail | grep Location | awk -F' ' '{print $2 "/hail"}')
export SPARK_HOME=/opt/conda/lib/python3.7/site-packages/pyspark
export JAVA_HOME=/opt/conda
time python /storage1/fs1/jin810/Active/fup/test_pyspark_02/test_hail_jointCalling.py
```

[[back to outline]](#outline)

Hail PYTHON SCRIPT: `/storage1/fs1/jin810/Active/fup/test_pyspark_02/test_hail_jointCalling.py`

```python
import os
import glob
import time
import pyspark
import hail as hl

start_time = time.time()

threads = int(os.environ['LSB_MAX_NUM_PROCESSORS']) - 1
hail_jars = "/opt/conda/lib/python3.7/site-packages/hail/backend/hail-all-spark.jar"
conf = pyspark.SparkConf().setAll([
    ('spark.master', 'local[{}]'.format(threads)),
    ('spark.app.name', 'Hail'),
    ('spark.jars', str(hail_jars)),
    ('spark.driver.extraClassPath', str(hail_jars)),
    ('spark.executor.extraClassPath', './hail-all-spark.jar'),
    ('spark.serializer', 'org.apache.spark.serializer.KryoSerializer'),
    ('spark.kryo.registrator', 'is.hail.kryo.HailKryoRegistrator'),
    ### https://discuss.hail.is/t/turning-run-combiner-performance-for-hail-local-mode/2318
    ('spark.driver.memory', '10g'),
    ('spark.executor.memory', '10g'),
    ])

### Using sc:
sc = pyspark.SparkContext(conf=conf)
hl.init(default_reference='GRCh38',sc=sc) 

### Input:
project_folder = "/storage1/fs1/jin810/Active/Neuropathy_WGS_2021May"
splitGVCF_Folder = project_folder + "/testHail_JointCalling/chr5_gvcf"
inputs = glob.glob(splitGVCF_Folder+"/*_germline.g.vcf.bgz")
### Output:
output_folder = "/storage1/fs1/jin810/Active/fup/test_pyspark_02"
temp_folder = output_folder + "/temp"
os.makedirs(temp_folder, mode=777, exist_ok=True)
output_file = output_folder + "/test_hail_run_combiner_chr5_v2.mt"

hl.experimental.run_combiner(inputs, 
                             use_genome_default_intervals=True,
                             out_file=output_file, 
                             tmp_path=temp_folder,
                             overwrite=True,
                             reference_genome='GRCh38')

processing_time_in_seconds = time.time() - start_time
print("--- %s minute ---" % (processing_time_in_seconds/60))
```

[[back to outline]](#outline)
