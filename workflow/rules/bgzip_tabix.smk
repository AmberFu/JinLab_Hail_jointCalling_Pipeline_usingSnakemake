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
