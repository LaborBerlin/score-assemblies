from glob import glob
shell.executable("/bin/bash")

snake_dir = workflow.basedir
report_rmd = snake_dir + "/scripts/report.Rmd"

workflow.global_resources['wget_busco'] = 1

busco_lineage = config.get('busco_lineage', 'bacteria')
run_bakta = config.get('bakta', '0')

wildcard_constraints:
  id = "[^/\\\\]+",
  ref = "[^/\\\\]+"

assemblies, = glob_wildcards("assemblies/{id,[^/\\\\]+}.fa")

if len(assemblies) == 0:
	print("Found no *.fa files in folder assemblies/")
	quit()

references, = glob_wildcards("references/{ref,[^/\\\\]+}.fa")
references_protein, = glob_wildcards("references-protein/{ref,[^/\\\\]+}.faa")

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
	list_assess_assembly_summ = expand("pomoxis/{id}/assess_assembly/{id}_{ref}_summ.txt", id=assemblies, ref=references)
	list_assess_assembly_meanQ_tsv = expand("pomoxis/{id}/assess_assembly/{id}_{ref}_scores.tsv", id=assemblies, ref=references)
	list_assess_assembly_meanQ_pdf = expand("pomoxis/{ref}_assess_assembly_all_meanQ.pdf", ref=references)
	list_assess_homopolymers_rel_len = expand("pomoxis/{id}/assess_homopolymers/{id}_{ref}_rel_len.tsv", id=assemblies, ref=references)
	list_assess_homopolymers_correct_len_tsv = expand("pomoxis/{id}/assess_homopolymers/{id}_{ref}_correct_len.tsv", id=assemblies, ref=references)
	list_assess_homopolymers_correct_len_pdf = expand("pomoxis/{ref}_assess_homopolymers_all_correct_len.pdf", ref=references)
	list_pomoxis = ["pomoxis/assess_assembly_all_scores.tsv", "pomoxis/assess_homopolymers_all_rel_len.tsv", "pomoxis/assess_homopolymers_all_correct_len.tsv" ]
	list_quast_report = expand("quast/{ref}/report.html", ref=references)
	list_dnadiff_report = expand("dnadiff/{ref}/{id}-dnadiff.report", id=assemblies, ref=references)
	list_dnadiff = [ "dnadiff/all_stats.tsv" ]
	list_dnadiff_tsv = expand("dnadiff/{ref}/{id}-dnadiff-stats.tsv", id=assemblies, ref=references)
	list_dnadiff_pdf = expand("dnadiff/{ref}_dnadiff_stats.pdf", ref=references)
	list_nucdiff_stat = expand("nucdiff/{ref}/{id}-nucdiff/results/nucdiff_stat.out", id=assemblies, ref=references)
	list_nucdiff_tsv = expand("nucdiff/{ref}/{id}-nucdiff/nucdiff.tsv", id=assemblies, ref=references)
	list_nucdiff = [ "nucdiff/all_stats.tsv" ]
	list_nucdiff_pdf = expand("nucdiff/{ref}_nucdiff_stats.pdf", ref=references)

if len(references_protein) > 0:
	list_ideel_ref_tsv = expand("ideel/diamond-ref/{ref}/{id}_{ref}.tsv", id=assemblies, ref=references_protein)
	list_ideel_ref_pdf = expand("ideel/{ref}_ideel_stats.pdf", ref=references_protein)

list_busco_out = expand("busco/{id}/short_summary.specific.{blin}_odb10.{id}.txt", id=assemblies, blin=busco_lineage)
list_busco_tsv = expand("busco/{id}/short_summary.specific.{blin}_odb10.{id}.tsv", id=assemblies, blin=busco_lineage)

list_prodigal_proteins = expand("ideel/prodigal/{id}.faa", id=assemblies)
list_ideel_uniprot_tsv = expand("ideel/diamond/{id}.tsv", id=assemblies)

list_bakta_out = []
if run_bakta == 1:
	list_bakta_out = expand("bakta/{id}/{id}.txt", id=assemblies)

rule all:
	input:
		list_assess_assembly_summ,
		list_assess_assembly_meanQ_tsv,
		list_assess_assembly_meanQ_pdf,
		list_assess_homopolymers_rel_len,
		list_assess_homopolymers_correct_len_tsv,
		list_assess_homopolymers_correct_len_pdf,
		list_pomoxis,
		list_busco_out,
		list_busco_tsv,
		"busco/all_stats.tsv",
		"busco/busco_stats.pdf",
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
		"ideel/uniprot/uniprot_sprot.fasta.gz",
		"ideel/uniprot/uniprot_sprot.dmnd",
		list_prodigal_proteins,
		list_ideel_uniprot_tsv,
		"ideel/ideel_stats.pdf",
		list_ideel_ref_tsv,
		list_ideel_ref_pdf,
		"report.html"

# -------------------------------------------------------------------------------------------------------------------------------------------
# pomoxis
# -------------------------------------------------------------------------------------------------------------------------------------------

rule assess_assembly:
	input:
		assembly = "assemblies/{id}.fa",
		ref = "references/{ref}.fa"
	output:
		summ = "pomoxis/{id}/assess_assembly/{id}_{ref}_summ.txt",
		tsv = "pomoxis/{id}/assess_assembly/{id}_{ref}_scores.tsv"
	log: "log/pomoxis/{id}/assess_assembly/{id}_{ref}_log.txt"
	shell:
		"""
		assess_assembly -r {input.ref} -i {input.assembly} -p pomoxis/{wildcards.id}/assess_assembly/{wildcards.id}_{wildcards.ref} >{log} 2>&1
		meanQ=$(grep -A2 '#  Q Scores' {output.summ} | tail -n1 | awk '{{print $2}}')
		percErr=$(grep -A2 '#  Percentage Errors' {output.summ} | tail -n1 | awk '{{print $2}}')
		echo "{wildcards.id}\t{wildcards.ref}\t$meanQ\t$percErr" > {output.tsv}
		"""

rule assess_homopolymers_minimap:
	threads: 5
	input:
		assembly = "assemblies/{id}.fa",
		ref = "references/{ref}.fa"
	output: "pomoxis/{id}/assess_homopolymers/{id}_{ref}.bam"
	log: "log/pomoxis/{id}/assess_homopolymers/{id}_{ref}.minimap2.log"
	shell:
		"""
		minimap2 -x asm5 -t {threads} --MD -a {input.ref} {input.assembly} 2>{log} | samtools sort -o {output} - >>{log} 2>&1
		"""

rule assess_homopolymers:
	input: "pomoxis/{id}/assess_homopolymers/{id}_{ref}.bam"
	output:
		count = "pomoxis/{id}/assess_homopolymers/{id}_{ref}_count/hp_counts.pkl",
		rel_len = "pomoxis/{id}/assess_homopolymers/{id}_{ref}_analyse/hp_rel_len_counts.txt",
		rel_len2 = "pomoxis/{id}/assess_homopolymers/{id}_{ref}_rel_len.tsv",
		correct_len = "pomoxis/{id}/assess_homopolymers/{id}_{ref}_analyse/hp_correct_vs_len.txt",
		correct_len2 = "pomoxis/{id}/assess_homopolymers/{id}_{ref}_correct_len.tsv"
	params:
		count_dir = directory("pomoxis/{id}/assess_homopolymers/{id}_{ref}_count"),
		analyse_dir = directory("pomoxis/{id}/assess_homopolymers/{id}_{ref}_analyse")
	log: "log/pomoxis/{id}/assess_homopolymers/{id}_{ref}_log.txt"
	shell:
		"""
		rm -rf {params.count_dir} {params.analyse_dir}
		assess_homopolymers count -o {params.count_dir} {input} >{log} 2>&1
		assess_homopolymers analyse -o {params.analyse_dir} {output.count} >>{log} 2>&1

		grep -v count {output.rel_len} | sed 's/^/{wildcards.id}\t{wildcards.ref}\t/' > {output.rel_len2}
		grep -v rlen {output.correct_len} | sed 's/^/{wildcards.id}\t{wildcards.ref}\t/' > {output.correct_len2}
		"""

rule gather_stats_pomoxis:
	input:
		aa_meanQ = list_assess_assembly_meanQ_tsv,
		hp_rel_len = list_assess_homopolymers_rel_len,
		hp_correct_len = list_assess_homopolymers_correct_len_tsv
	output:
		all_aa_meanQ = "pomoxis/assess_assembly_all_scores.tsv",
		all_hp_rel_len = "pomoxis/assess_homopolymers_all_rel_len.tsv",
		all_hp_correct_len = "pomoxis/assess_homopolymers_all_correct_len.tsv"
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
	threads: 5
	resources: wget_busco=1
	input:
		assembly = "assemblies/{id}.fa"
	output:
		"busco/{id}/short_summary.specific.{busco_lineage}_odb10.{id}.txt",
	log: "log/busco/{id}/busco_{busco_lineage}.log"
	shell:
		"""
		cd busco && busco -q -c {threads} -f -m genome -l {busco_lineage} -o {wildcards.id} -i ../{input} >../{log} 2>&1
		"""

rule busco2tsv:
	threads: 1
	input:
		"busco/{id}/short_summary.specific.{busco_lineage}_odb10.{id}.txt",
	output:
		"busco/{id}/short_summary.specific.{busco_lineage}_odb10.{id}.tsv"
	shell:
		"""
		perl -lne 'print "{wildcards.id}\t$1\t$2\t$3\t$4" if /C:([\-\d\.]+).*F:([\-\d\.]+).*M:([\-\d\.]+).*n:(\d+)/' {input} > {output}
		"""

rule gather_stats_busco:
	input:
		list_busco_tsv
	output:
		"busco/all_stats.tsv"
	shell:
		"""
		cat {input} > {output}
		"""

# -------------------------------------------------------------------------------------------------------------------------------------------
# QUAST
# -------------------------------------------------------------------------------------------------------------------------------------------

rule quast:
	threads: 2
	input:
		reference = "references/{ref}.fa"
	output:
		report = "quast/{ref}/report.html"
	log: "log/quast/{ref}/quast.log"
	shell:
		"""
		quast -t {threads} --glimmer -o quast/{wildcards.ref} -r {input.reference} assemblies/*.fa >{log} 2>&1
		"""

# -------------------------------------------------------------------------------------------------------------------------------------------
# dnadiff
# -------------------------------------------------------------------------------------------------------------------------------------------

rule dnadiff:
	threads: 1
	input:
		reference = "references/{ref}.fa",
		assembly = "assemblies/{id}.fa"
	output:
		dnadiff_report = "dnadiff/{ref}/{id}-dnadiff.report",
		stats_tsv = "dnadiff/{ref}/{id}-dnadiff-stats.tsv"
	log: "log/dnadiff/{ref}/{id}-dnadiff.log"
	shell:
		"""
		dnadiff -p dnadiff/{wildcards.ref}/{wildcards.id}-dnadiff {input.reference} {input.assembly} >{log} 2>&1
		cat {output.dnadiff_report} | grep -A3 '1-to-1' | grep 'AvgIdentity' | sed -e 's/^/{wildcards.id}\t{wildcards.ref}\t/' | perl -lpne 's/\s+/\t/g' > {output.stats_tsv}
		grep TotalIndels {output.dnadiff_report} | sed -e 's/^/{wildcards.id}\t{wildcards.ref}\t/' | perl -lpne 's/\s+/\t/g' >> {output.stats_tsv}
		"""

rule gather_stats_dnadiff:
	input:
		all_reports = list_dnadiff_tsv
	output:
		all_dnadiff_report_tsv = "dnadiff/all_stats.tsv"
	shell:
		"""
		cat {input} > {output}
		"""

# -------------------------------------------------------------------------------------------------------------------------------------------
# nucdiff
# -------------------------------------------------------------------------------------------------------------------------------------------

rule nucdiff:
	threads: 1
	input:
		reference = "references/{ref}.fa",
		assembly = "assemblies/{id}.fa"
	output:
		nucdiff_stat = "nucdiff/{ref}/{id}-nucdiff/results/nucdiff_stat.out",
		nucdiff_tsv = "nucdiff/{ref}/{id}-nucdiff/nucdiff.tsv"
	params:
		nucdiff_dir = directory("nucdiff/{ref}/{id}-nucdiff/")
	log: "log/nucdiff/{ref}/{id}-nucdiff.log"
	shell:
		"""
		nucdiff {input.reference} {input.assembly} {params.nucdiff_dir} nucdiff >{log} 2>&1
		cat {output.nucdiff_stat} | grep 'Insertions\|Deletions\|Substitutions' | sed -e 's/^/{wildcards.id}\t{wildcards.ref}\t/' > {output.nucdiff_tsv}
		"""

rule gather_stats_nucdiff:
	input:
		all_reports = list_nucdiff_tsv
	output:
		all_nucdiff_tsv = "nucdiff/all_stats.tsv"
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
		"ideel/uniprot/uniprot_sprot.fasta.gz"
	log: "log/ideel/uniprot/download.log"
	shell:
		"""
		wget -P ideel/uniprot ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz >{log} 2>&1
		"""

rule diamond_makedb:
	input:
		"ideel/uniprot/uniprot_sprot.fasta.gz"
	output:
		"ideel/uniprot/uniprot_sprot.dmnd"
	log: "log/diamond-makedb-uniprot.log"
	shell:
		"""
		diamond makedb --db ideel/uniprot/uniprot_sprot --in ideel/uniprot/uniprot_sprot.fasta.gz >{log} 2>&1
		"""

rule prodigal:
	input: "assemblies/{id}.fa"
	output: "ideel/prodigal/{id}.faa"
	log: "log/ideel/prodigal/{id}.log"
	shell:
		"""
		prodigal -a {output} -i {input} >{log} 2>&1
		# remove stop codon (*) from end of sequences
		sed -i 's/*$//' {output}
		"""

rule diamond:
  threads: 5
	input:
		proteins = "ideel/prodigal/{id}.faa",
		db = "ideel/uniprot/uniprot_sprot.dmnd"
	output: "ideel/diamond/{id}.tsv"
	log: "log/ideel/diamond/{id}.log"
	shell:
		"""
		diamond blastp --threads {threads} --max-target-seqs 1 --db {input.db} --query {input.proteins} --outfmt 6 qlen slen --out {output} >{log} 2>&1
		"""

# -------------------------------------------------------------------------------------------------------------------------------------------
# bakta
# -------------------------------------------------------------------------------------------------------------------------------------------
rule download_bakta_db:
	conda: "env/env-bakta.yaml"
	priority: 10
	output: "bakta/db/version.json"
	log: "log/bakta/download.log"
	shell:
		"""
		wget -N -P bakta https://jlubox.uni-giessen.de/dl/fiWhS6LJi8AizXvspaRxRzPN/db.tar.gz >{log} 2>&1
		tar --directory bakta -xf bakta/db.tar.gz >{log} 2>&1
		amrfinder_update --force_update --database bakta/db/amrfinderplus-db/ >{log} 2>&1
		"""

rule bakta:
	conda: "env/env-bakta.yaml"
	threads: 5
	input:
		fa = "assemblies/{id}.fa",
		db = "bakta/db/version.json"
	output: "bakta/{id}/{id}.txt"
	log: "log/bakta/{id}.log"
	shell:
		"""
		bakta --db bakta/db --verbose --output bakta/{wildcards.id} --prefix {wildcards.id} --threads {threads} {input.fa} >{log} 2>&1
		"""
# -------------------------------------------------------------------------------------------------------------------------------------------
# protein comparison to reference genome
# -------------------------------------------------------------------------------------------------------------------------------------------

rule diamond_ref_makedb:
	input:
		"references-protein/{ref}.faa"
	output:
		"references-protein/{ref}.dmnd"
	log: "log/references-protein/{ref}-diamond-makedb.log"
	shell:
		"""
		diamond makedb --db {output} --in {input} >{log} 2>&1
		"""

rule diamond_ref:
  threads: 5
	input:
		proteins = "ideel/prodigal/{id}.faa",
		db = "references-protein/{ref}.dmnd"
	output: "ideel/diamond-ref/{ref}/{id}_{ref}.tsv"
	log: "log/ideel/diamond-ref/{ref}_{id}.log"
	shell:
		"""
		diamond blastp --threads {threads} --max-target-seqs 1 --db {input.db} --query {input.proteins} --outfmt 6 qlen slen --out {output} >{log} 2>&1
		"""
# -------------------------------------------------------------------------------------------------------------------------------------------
# plots
# -------------------------------------------------------------------------------------------------------------------------------------------

rule plot_assess_assembly:
	input:
		"pomoxis/assess_assembly_all_scores.tsv"
	output:
		list_assess_assembly_meanQ_pdf
	script: "scripts/plot-assess_assembly.R"

rule plot_asses_homopolymers:
	input:
		"pomoxis/assess_homopolymers_all_correct_len.tsv"
	output:
		list_assess_homopolymers_correct_len_pdf,
	script: "scripts/plot-assess_homopolymers.R"

rule plot_busco:
	input:
		"busco/all_stats.tsv"
	output:
		"busco/busco_stats.pdf",
	script: "scripts/plot-busco.R"

rule plot_dnadiff:
	input:
		"dnadiff/all_stats.tsv"
	output:
		list_dnadiff_pdf,
	script: "scripts/plot-dnadiff.R"

rule plot_nucdiff:
	input:
		"nucdiff/all_stats.tsv"
	output:
		list_nucdiff_pdf
	script: "scripts/plot-nucdiff.R"

rule plot_ideel:
	input:
		list_ideel_ref_tsv,
		list_ideel_uniprot_tsv
	output:
		"ideel/ideel_stats.pdf",
		list_ideel_ref_pdf
	script: "scripts/plot-ideel.R"

rule report_html:
	input:
		"busco/all_stats.tsv",
		list_ideel_uniprot_tsv,
		list_ideel_ref_tsv,
		list_pomoxis
	output:
		"report.html"
	message:
		"knit report HTML"
	log:
		"log/report.log"
	params:
		wd = os.getcwd()
	shell:
		"""
		Rscript -e 'args<-commandArgs(trailingOnly = TRUE); rmarkdown::render(args[1], output_file=args[3], knit_root_dir=args[4])' {report_rmd} {snake_dir} {params.wd}/{output} {params.wd} >{log} 2>&1
		"""

