from glob import glob


out_dir = "score-assemblies-data"
log_dir = "score-assemblies-data/log"

report_rmd = workflow.source_path("scripts/report.Rmd")

workflow.global_resources["wget_busco"] = 1

busco_lineage = config.get("busco_lineage", "bacteria")
run_bakta = config.get("bakta", "0")


wildcard_constraints:
    id="[^/\\\\]+",
    ref="[^/\\\\]+",


(assemblies,) = glob_wildcards("assemblies/{id,[^/\\\\]+}.fa")

assemblies_fa = expand("assemblies/{id}.fa", id=assemblies)

if len(assemblies) == 0:
    print("Found no *.fa files in folder assemblies/")
    quit()

(references,) = glob_wildcards("references/{ref,[^/\\\\]+}.fa")
(references_protein,) = glob_wildcards("references-protein/{ref,[^/\\\\]+}.faa")

list_assess_assembly_summ = []
list_assess_assembly_meanQ_tsv = []
list_assess_assembly_meanQ_pdf = []
list_assess_homopolymers_rel_len = []
list_assess_homopolymers_correct_len_tsv = []
list_assess_homopolymers_correct_len_pdf = []
list_pomoxis = []
list_quast_report = []
list_dnadiff_report = []
list_dnadiff = []
list_dnadiff_tsv = []
list_dnadiff_pdf = []
list_nucdiff_stat = []
list_nucdiff_tsv = []
list_nucdiff = []
list_nucdiff_pdf = []
list_ideel_ref_tsv = []
list_ideel_ref_pdf = []
if len(references) > 0:
    list_assess_assembly_summ = expand(
        out_dir + "/pomoxis/{id}/assess_assembly/{id}_{ref}_summ.txt",
        id=assemblies,
        ref=references,
    )
    list_assess_assembly_meanQ_tsv = expand(
        out_dir + "/pomoxis/{id}/assess_assembly/{id}_{ref}_scores.tsv",
        id=assemblies,
        ref=references,
    )
    list_assess_assembly_meanQ_pdf = expand(
        out_dir + "/pomoxis/{ref}_assess_assembly_all_meanQ.pdf", ref=references
    )
    list_assess_homopolymers_rel_len = expand(
        out_dir + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_rel_len.tsv",
        id=assemblies,
        ref=references,
    )
    list_assess_homopolymers_correct_len_tsv = expand(
        out_dir + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_correct_len.tsv",
        id=assemblies,
        ref=references,
    )
    list_pomoxis = [
        out_dir + "/pomoxis/assess_assembly_all_scores.tsv",
        out_dir + "/pomoxis/assess_homopolymers_all_rel_len.tsv",
        out_dir + "/pomoxis/assess_homopolymers_all_correct_len.tsv",
    ]
    list_quast_report = expand(out_dir + "/quast/{ref}/report.html", ref=references)
    list_dnadiff_report = expand(
        out_dir + "/dnadiff/{ref}/{id}-dnadiff.report", id=assemblies, ref=references
    )
    list_dnadiff = [out_dir + "/dnadiff/all_stats.tsv"]
    list_dnadiff_tsv = expand(
        out_dir + "/dnadiff/{ref}/{id}-dnadiff-stats.tsv", id=assemblies, ref=references
    )
    list_dnadiff_pdf = expand(
        out_dir + "/dnadiff/{ref}_dnadiff_stats.pdf", ref=references
    )
    list_nucdiff_stat = expand(
        out_dir + "/nucdiff/{ref}/{id}-nucdiff/results/nucdiff_stat.out",
        id=assemblies,
        ref=references,
    )
    list_nucdiff_tsv = expand(
        out_dir + "/nucdiff/{ref}/{id}-nucdiff/nucdiff.tsv",
        id=assemblies,
        ref=references,
    )
    list_nucdiff = [out_dir + "/nucdiff/all_stats.tsv"]
    list_nucdiff_pdf = expand(
        out_dir + "/nucdiff/{ref}_nucdiff_stats.pdf", ref=references
    )

if len(references_protein) > 0:
    list_ideel_ref_tsv = expand(
        out_dir + "/ideel/diamond-ref/{ref}/{id}_{ref}.tsv",
        id=assemblies,
        ref=references_protein,
    )
    list_ideel_ref_pdf = expand(
        out_dir + "/ideel/{ref}_ideel_{type}.pdf",
        ref=references_protein,
        type=["histograms", "boxplots"],
    )

list_busco_out = expand(
    out_dir + "/busco/{id}/short_summary.specific.{blin}_odb10.{id}.txt",
    id=assemblies,
    blin=busco_lineage,
)
list_busco_tsv = expand(
    out_dir + "/busco/{id}/short_summary.specific.{blin}_odb10.{id}.tsv",
    id=assemblies,
    blin=busco_lineage,
)

list_prodigal_proteins = expand(out_dir + "/ideel/prodigal/{id}.faa", id=assemblies)
list_ideel_uniprot_tsv = expand(out_dir + "/ideel/diamond/{id}.tsv", id=assemblies)

list_bakta_out = []
if run_bakta == 1:
    list_bakta_out = expand(out_dir + "/bakta/{id}/{id}.txt", id=assemblies)


rule all:
    input:
        list_assess_assembly_summ,
        list_assess_assembly_meanQ_tsv,
        list_assess_assembly_meanQ_pdf,
        list_assess_homopolymers_rel_len,
        list_assess_homopolymers_correct_len_tsv,
        list_pomoxis,
        list_busco_out,
        list_busco_tsv,
        out_dir + "/busco/all_stats.tsv",
        out_dir + "/busco/busco_stats.pdf",
        list_bakta_out,
        list_quast_report,
        list_dnadiff_report,
        list_dnadiff,
        list_dnadiff_tsv,
        list_dnadiff_pdf,
        list_nucdiff_stat,
        list_nucdiff_tsv,
        list_nucdiff,
        list_nucdiff_pdf,
        out_dir + "/ideel/uniprot/uniprot_sprot.fasta.gz",
        out_dir + "/ideel/uniprot/uniprot_sprot.dmnd",
        list_prodigal_proteins,
        out_dir + "/ideel/prodigal_stats.tsv",
        list_ideel_uniprot_tsv,
        out_dir + "/ideel/ideel_uniprot_histograms.pdf",
        out_dir + "/ideel/ideel_uniprot_boxplots.pdf",
        list_ideel_ref_tsv,
        list_ideel_ref_pdf,
        "score-assemblies-report.html",


# -------------------------------------------------------------------------------------------------------------------------------------------
# pomoxis
# -------------------------------------------------------------------------------------------------------------------------------------------


rule assess_assembly:
    conda:
        "env/env-pomoxis.yaml"
    input:
        assembly="assemblies/{id}.fa",
        ref="references/{ref}.fa",
    output:
        summ=out_dir + "/pomoxis/{id}/assess_assembly/{id}_{ref}_summ.txt",
        tsv=out_dir + "/pomoxis/{id}/assess_assembly/{id}_{ref}_scores.tsv",
    log:
        log_dir + "/pomoxis/{id}/assess_assembly/{id}_{ref}_log.txt",
    params:
        out_dir=out_dir,
    shell:
        """
        assess_assembly -r {input.ref} -i {input.assembly} -p {params.out_dir}/pomoxis/{wildcards.id}/assess_assembly/{wildcards.id}_{wildcards.ref} >{log} 2>&1
        meanQ=$(grep -A2 '#  Q Scores' {output.summ} | tail -n1 | awk '{{print $2}}')
        percErr=$(grep -A2 '#  Percentage Errors' {output.summ} | tail -n1 | awk '{{print $2}}')
        echo "{wildcards.id}\t{wildcards.ref}\t$meanQ\t$percErr" > {output.tsv}
        """


rule assess_homopolymers_minimap:
    conda:
        "env/env-pomoxis.yaml"
    threads: 5
    input:
        assembly="assemblies/{id}.fa",
        ref="references/{ref}.fa",
    output:
        out_dir + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}.bam",
    log:
        log_dir + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}.minimap2.log",
    shell:
        """
        minimap2 -x asm5 -t {threads} --MD -a {input.ref} {input.assembly} 2>{log} | samtools sort -o {output} - >>{log} 2>&1
        """


rule assess_homopolymers:
    conda:
        "env/env-pomoxis.yaml"
    input:
        out_dir + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}.bam",
    output:
        hp_count=out_dir
        + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_count/hp_counts.pkl",
        rel_len=out_dir
        + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_analyse/hp_rel_len_counts.txt",
        rel_len2=out_dir + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_rel_len.tsv",
        correct_len=out_dir
        + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_analyse/hp_correct_vs_len.txt",
        correct_len2=out_dir
        + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_correct_len.tsv",
    params:
        count_dir=directory(
            out_dir + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_count"
        ),
        analyse_dir=directory(
            out_dir + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_analyse"
        ),
    log:
        log_dir + "/pomoxis/{id}/assess_homopolymers/{id}_{ref}_log.txt",
    shell:
        """
        rm -rf {params.count_dir} {params.analyse_dir}
        assess_homopolymers count -o {params.count_dir} {input} >{log} 2>&1
        assess_homopolymers analyse -o {params.analyse_dir} {output.hp_count} >>{log} 2>&1

        grep -v count {output.rel_len} | sed 's/^/{wildcards.id}\t{wildcards.ref}\t/' > {output.rel_len2}
        grep -v rlen {output.correct_len} | sed 's/^/{wildcards.id}\t{wildcards.ref}\t/' > {output.correct_len2}
        """


rule gather_stats_pomoxis:
    input:
        aa_meanQ=list_assess_assembly_meanQ_tsv,
        hp_rel_len=list_assess_homopolymers_rel_len,
        hp_correct_len=list_assess_homopolymers_correct_len_tsv,
    output:
        all_aa_meanQ=out_dir + "/pomoxis/assess_assembly_all_scores.tsv",
        all_hp_rel_len=out_dir + "/pomoxis/assess_homopolymers_all_rel_len.tsv",
        all_hp_correct_len=out_dir + "/pomoxis/assess_homopolymers_all_correct_len.tsv",
    shell:
        """
        #filter inf values, which happen when comparing identical assemblies
        cat {input.aa_meanQ} | grep -v -w 'inf' > {output.all_aa_meanQ}

        cat {input.hp_rel_len}  > {output.all_hp_rel_len}
        cat {input.hp_correct_len}  > {output.all_hp_correct_len}
        """


# -------------------------------------------------------------------------------------------------------------------------------------------
# BUSCO
# -------------------------------------------------------------------------------------------------------------------------------------------


rule busco:
    conda:
        "env/env-busco.yaml"
    threads: 5
    resources:
        wget_busco=1,
    input:
        assembly="assemblies/{id}.fa",
    output:
        out_dir + "/busco/{id}/short_summary.specific." + busco_lineage + "_odb10.{id}.txt",
    params:
        out_dir=out_dir,
        busco_lineage=busco_lineage,
    log:
        log_dir + "/busco/{id}/busco_" + busco_lineage + ".log",
    shell:
        """
        cd {params.out_dir}/busco && busco -q -c {threads} -f -m genome -l {params.busco_lineage} -o {wildcards.id} -i ../../{input} >../../{log} 2>&1
        """


rule busco2tsv:
    threads: 1
    input:
        out_dir + "/busco/{id}/short_summary.specific.{busco_lineage}_odb10.{id}.txt",
    output:
        out_dir + "/busco/{id}/short_summary.specific.{busco_lineage}_odb10.{id}.tsv",
    shell:
        """
        perl -lne 'print "{wildcards.id}\t$1\t$2\t$3\t$4" if /C:([\-\d\.]+).*F:([\-\d\.]+).*M:([\-\d\.]+).*n:(\d+)/' {input} > {output}
        """


rule gather_stats_busco:
    input:
        list_busco_tsv,
    output:
        out_dir + "/busco/all_stats.tsv",
    shell:
        """
        cat {input} > {output}
        """


# -------------------------------------------------------------------------------------------------------------------------------------------
# QUAST
# -------------------------------------------------------------------------------------------------------------------------------------------


rule quast:
    conda:
        "env/env-quast.yaml"
    threads: 5
    input:
        reference="references/{ref}.fa",
        fa=assemblies_fa,
    output:
        report=out_dir + "/quast/{ref}/report.html",
    log:
        log_dir + "/quast/{ref}/quast.log",
    params:
        out_dir=out_dir,
    shell:
        """
        quast -t {threads} --glimmer -o {params.out_dir}/quast/{wildcards.ref} -r {input.reference} {input.fa} >{log} 2>&1
        """


# -------------------------------------------------------------------------------------------------------------------------------------------
# dnadiff
# -------------------------------------------------------------------------------------------------------------------------------------------


rule dnadiff:
    conda:
        "env/env-mummer.yaml"
    threads: 1
    input:
        reference="references/{ref}.fa",
        assembly="assemblies/{id}.fa",
    output:
        dnadiff_report=out_dir + "/dnadiff/{ref}/{id}-dnadiff.report",
        stats_tsv=out_dir + "/dnadiff/{ref}/{id}-dnadiff-stats.tsv",
    log:
        log_dir + "/dnadiff/{ref}/{id}-dnadiff.log",
    params:
        out_dir=out_dir,
    shell:
        """
        perl $CONDA_PREFIX/bin/dnadiff -p {params.out_dir}/dnadiff/{wildcards.ref}/{wildcards.id}-dnadiff {input.reference} {input.assembly} >{log} 2>&1
        cat {output.dnadiff_report} | grep -A3 '1-to-1' | grep 'AvgIdentity' | sed -e 's/^/{wildcards.id}\t{wildcards.ref}\t/' | perl -lpne 's/\s+/\t/g' > {output.stats_tsv}
        grep TotalIndels {output.dnadiff_report} | sed -e 's/^/{wildcards.id}\t{wildcards.ref}\t/' | perl -lpne 's/\s+/\t/g' >> {output.stats_tsv}
        """


rule gather_stats_dnadiff:
    input:
        all_reports=list_dnadiff_tsv,
    output:
        all_dnadiff_report_tsv=out_dir + "/dnadiff/all_stats.tsv",
    shell:
        """
        cat {input} > {output}
        """


# -------------------------------------------------------------------------------------------------------------------------------------------
# nucdiff
# -------------------------------------------------------------------------------------------------------------------------------------------


rule nucdiff:
    conda:
        "env/env-nucdiff.yaml"
    threads: 1
    input:
        reference="references/{ref}.fa",
        assembly="assemblies/{id}.fa",
    output:
        nucdiff_stat=out_dir + "/nucdiff/{ref}/{id}-nucdiff/results/nucdiff_stat.out",
        nucdiff_tsv=out_dir + "/nucdiff/{ref}/{id}-nucdiff/nucdiff.tsv",
    params:
        nucdiff_dir=directory(out_dir + "/nucdiff/{ref}/{id}-nucdiff/"),
    log:
        log_dir + "/nucdiff/{ref}/{id}-nucdiff.log",
    shell:
        """
        nucdiff {input.reference} {input.assembly} {params.nucdiff_dir} nucdiff >{log} 2>&1
        cat {output.nucdiff_stat} | grep 'Insertions\|Deletions\|Substitutions' | sed -e 's/^/{wildcards.id}\t{wildcards.ref}\t/' > {output.nucdiff_tsv}
        """


rule gather_stats_nucdiff:
    input:
        all_reports=list_nucdiff_tsv,
    output:
        all_nucdiff_tsv=out_dir + "/nucdiff/all_stats.tsv",
    shell:
        """
        cat {input} > {output}
        """


# -------------------------------------------------------------------------------------------------------------------------------------------
# ideel with uniprot
# -------------------------------------------------------------------------------------------------------------------------------------------


rule download_uniprot:
    priority: 10
    output:
        out_dir + "/ideel/uniprot/uniprot_sprot.fasta.gz",
    log:
        log_dir + "/ideel/uniprot/download.log",
    params:
      out_dir = out_dir
    shell:
        """
        wget -P {params.out_dir}/ideel/uniprot http://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz >{log} 2>&1
        """


rule diamond_makedb:
    conda:
        "env/env-ideel.yaml"
    input:
        out_dir + "/ideel/uniprot/uniprot_sprot.fasta.gz",
    output:
        out_dir + "/ideel/uniprot/uniprot_sprot.dmnd",
    params:
        out_dir=out_dir,
    log:
        log_dir + "/diamond-makedb-uniprot.log",
    shell:
        """
        diamond makedb --db {params.out_dir}/ideel/uniprot/uniprot_sprot --in {input} >{log} 2>&1
        """


rule prodigal:
    conda:
        "env/env-ideel.yaml"
    input:
        "assemblies/{id}.fa",
    output:
        out_dir + "/ideel/prodigal/{id}.faa",
    log:
        log_dir + "/ideel/prodigal/{id}.log",
    shell:
        """
        prodigal -a {output} -i {input} >{log} 2>&1
        # remove stop codon (*) from end of sequences
        sed -i 's/*$//' {output}
        """


rule prodigal_stats:
    conda:
        "env/env-ideel.yaml"
    threads: 1
    input:
        list_prodigal_proteins,
    output:
        out_dir + "/ideel/prodigal_stats.tsv",
    log:
        log_dir + "/ideel/prodigal_stats.log",
    shell:
        """
        seqkit stats -T -t protein {input} >{output} 2>{log}
        """


rule diamond:
    conda:
        "env/env-ideel.yaml"
    threads: 5
    input:
        proteins=out_dir + "/ideel/prodigal/{id}.faa",
        db=out_dir + "/ideel/uniprot/uniprot_sprot.dmnd",
    output:
        out_dir + "/ideel/diamond/{id}.tsv",
    log:
        log_dir + "/ideel/diamond/{id}.log",
    shell:
        """
        diamond blastp --threads {threads} --max-target-seqs 1 --db {input.db} --query {input.proteins} --outfmt 6 qseqid sseqid qlen slen --out {output} >{log} 2>&1
        """


# -------------------------------------------------------------------------------------------------------------------------------------------
# bakta
# -------------------------------------------------------------------------------------------------------------------------------------------
rule download_bakta_db:
    conda:
        "env/env-bakta.yaml"
    priority: 10
    output:
        out_dir + "/bakta/db-light/version.json",
    params:
        out_dir=out_dir,
    log:
        log_dir + "/bakta/download.log",
    shell:
        """
        wget -N -P {params.out_dir}/bakta https://zenodo.org/record/7669534/files/db-light.tar.gz >{log} 2>&1
        tar --directory {params.out_dir}/bakta -xf {params.out_dir}/bakta/db-light.tar.gz >{log} 2>&1
        amrfinder_update --force_update --database {params.out_dir}/bakta/db-light/amrfinderplus-db/ >{log} 2>&1
        """


rule bakta:
    conda:
        "env/env-bakta.yaml"
    threads: 5
    input:
        fa="assemblies/{id}.fa",
        db=out_dir + "/bakta/db-light/version.json",
    output:
        out_dir + "/bakta/{id}/{id}.txt",
    params:
        out_dir=out_dir,
    log:
        log_dir + "/bakta/{id}.log",
    shell:
        """
        bakta --db {params.out_dir}/bakta/db-light --verbose --output {params.out_dir}/bakta/{wildcards.id} --prefix {wildcards.id} --threads {threads} {input.fa} >{log} 2>&1
        """


# -------------------------------------------------------------------------------------------------------------------------------------------
# protein comparison to reference genome
# -------------------------------------------------------------------------------------------------------------------------------------------


rule diamond_ref_makedb:
    conda:
        "env/env-ideel.yaml"
    input:
        "references-protein/{ref}.faa",
    output:
        "references-protein/{ref}.dmnd",
    log:
        log_dir + "/references-protein/{ref}-diamond-makedb.log",
    shell:
        """
        diamond makedb --db {output} --in {input} >{log} 2>&1
        """


rule diamond_ref:
    conda:
        "env/env-ideel.yaml"
    threads: 5
    input:
        proteins=out_dir + "/ideel/prodigal/{id}.faa",
        db="references-protein/{ref}.dmnd",
    output:
        out_dir + "/ideel/diamond-ref/{ref}/{id}_{ref}.tsv",
    log:
        log_dir + "/ideel/diamond-ref/{ref}_{id}.log",
    shell:
        """
        diamond blastp --threads {threads} --max-target-seqs 1 --db {input.db} --query {input.proteins} --outfmt 6 qseqid sseqid qlen slen --out {output} >{log} 2>&1
        """


# -------------------------------------------------------------------------------------------------------------------------------------------
# plots
# -------------------------------------------------------------------------------------------------------------------------------------------


rule plot_assess_assembly:
    conda:
        "env/env-r.yaml"
    input:
        out_dir + "/pomoxis/assess_assembly_all_scores.tsv",
    output:
        list_assess_assembly_meanQ_pdf,
    params:
        out_dir=out_dir,
    script:
        "scripts/plot-assess_assembly.R"


# rule plot_asses_homopolymers:
# 	conda: "env/env-r.yaml"
# 	input:
# 		"pomoxis/assess_homopolymers_all_correct_len.tsv"
# 	output:
# 		list_assess_homopolymers_correct_len_pdf,
# 	script: "scripts/plot-assess_homopolymers.R"


rule plot_busco:
    conda:
        "env/env-r.yaml"
    input:
        out_dir + "/busco/all_stats.tsv",
    output:
        out_dir + "/busco/busco_stats.pdf",
    params:
        out_dir=out_dir,
    script:
        "scripts/plot-busco.R"


rule plot_dnadiff:
    conda:
        "env/env-r.yaml"
    input:
        out_dir + "/dnadiff/all_stats.tsv",
    output:
        list_dnadiff_pdf,
    params:
        out_dir=out_dir,
    script:
        "scripts/plot-dnadiff.R"


rule plot_nucdiff:
    conda:
        "env/env-r.yaml"
    input:
        out_dir + "/nucdiff/all_stats.tsv",
    output:
        list_nucdiff_pdf,
    params:
        out_dir=out_dir,
    script:
        "scripts/plot-nucdiff.R"


rule plot_ideel:
    conda:
        "env/env-r.yaml"
    input:
        list_ideel_ref_tsv,
        list_ideel_uniprot_tsv,
    output:
        out_dir + "/ideel/ideel_uniprot_histograms.pdf",
        out_dir + "/ideel/ideel_uniprot_boxplots.pdf",
        list_ideel_ref_pdf,
    params:
        out_dir=out_dir,
    script:
        "scripts/plot-ideel.R"


rule report_html:
    conda:
        "env/env-r.yaml"
    input:
        out_dir + "/busco/all_stats.tsv",
        list_ideel_uniprot_tsv,
        list_ideel_ref_tsv,
        list_dnadiff,
        list_nucdiff,
        list_pomoxis,
        list_bakta_out,
    output:
        "score-assemblies-report.html",
    message:
        "report HTML"
    log:
        log_dir + "/report.log",
    params:
        wd=os.getcwd(),
        out_dir=out_dir,
        report_rmd=report_rmd,
    shell:
        """
        Rscript -e 'args<-commandArgs(trailingOnly = TRUE); rmarkdown::render(args[1], output_file=args[3], knit_root_dir=args[4])' {params.report_rmd} {params.out_dir} {params.wd}/{output} {params.wd} >{log} 2>&1
        """
