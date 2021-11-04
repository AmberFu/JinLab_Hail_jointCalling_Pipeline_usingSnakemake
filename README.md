# JinLab_Hail_jointCalling_Pipeline_usingSnakemake
TEST Snakemake Pipeline for Hail Joint-Calling (run_combiner()) Pipeline

Follow [snakemake - Best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html)

### Cluster setting

```json
{
    "__default__" :
    {
        "core" : 2,
        "mem" : "10GB",
        "resources" : "\"rusage[mem=10GB] span[hosts=1]\"",
        "image" : "spashleyfu/ubuntu20_snakemake:bamMetrics",
        "gpu" : "\"\"",
        "log" : "logs/out.txt",
        "err" : "logs/err.txt",
        "jobgroup": "/fup/jobGroup_50_snakemake"
    },
    "gvcf2gvcfbgz" : {
        "image" : "spashleyfu/ubuntu18_vep104:hail_gsutil",
        "log" : "logs/gvcf2gvcfbgz_out.txt",
        "err" : "logs/gvcf2gvcfbgz_err.txt"
    },
    "generate_tabix" : {
        "image" : "spashleyfu/ubuntu18_vep104:hail_gsutil",
        "log" : "logs/generate_tabix_out.txt",
        "err" : "logs/generate_tabix_err.txt"
    },
    "hail_run_combiner" : {
        "image" : "spashleyfu/ubuntu18_vep104:hail_gsutil",
        "log" : "logs/hail_run_combiner_out.txt",
        "err" : "logs/hail_run_combiner_err.txt"
    }
}
```


### Detail:

#### How to BGZ a GVCF file and TABIX it:

```
(base) fup@compute1-exec-283:/storage1/fs1/jin810/Active/Neuropathy_WGS_2021May$ time bgzip TWHJ-PNRR-10001-10001_1_germline.g.vcf && tabix -p vcf TWHJ-PNRR-10001-10001_1_germline.g.vcf.gz
real	1m40.566s
user	1m35.444s
sys	0m4.344s
```

#### How to run Hail Joint-Calling:

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
