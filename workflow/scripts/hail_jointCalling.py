import os
import time
import pyspark
import hail as hl

start_time = time.time()

# Set JAVA HOME
os.environ['JAVA_HOME'] = "/opt/conda"
os.environ['SPARK_HOME'] = "/opt/conda/lib/python3.7/site-packages/pyspark"
os.environ['HAIL_HOME'] = "/opt/conda/lib/python3.7/site-packages/hail"

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
