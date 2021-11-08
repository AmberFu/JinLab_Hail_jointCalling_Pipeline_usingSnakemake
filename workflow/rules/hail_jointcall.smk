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

