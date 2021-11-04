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
