shell.executable("/bin/bash")

busco_lineage = config.get('busco_lineage', 'bacteria')

assemblies, = glob_wildcards("assemblies/{id}.fa")
references, = glob_wildcards("references/{id}.fa")

list_assess_assembly_summ = expand("{id}/assess_assembly/{id}_{ref}_summ.txt", id=assemblies, ref=references)
list_assess_assembly_meanQ = expand("{id}/assess_assembly/{id}_{ref}_meanQ.tsv", id=assemblies, ref=references)
list_assess_homopolymers_rel_len = expand("{id}/assess_homopolymers/{id}_{ref}_rel_len.tsv", id=assemblies, ref=references)
list_assess_homopolymers_correct_len = expand("{id}/assess_homopolymers/{id}_{ref}_correct_len.tsv", id=assemblies, ref=references)
list_assess_homopolymers_correct_len_pdf = expand("plots/{ref}_assess_homopolymers_all_correct_len.pdf", ref=references)

list_busco_out = expand("{id}/busco/short_summary.specific.{blin}_odb10.busco.txt", id=assemblies, blin=busco_lineage)

rule all:
	input:
		list_assess_assembly_summ,
		list_assess_assembly_meanQ,
		expand("{id}/assess_homopolymers/{id}_{ref}.bam", id=assemblies, ref=references),
		expand("{id}/assess_homopolymers/{id}_{ref}_analyse/hp_rel_len_counts.txt", id=assemblies, ref=references),
		list_assess_homopolymers_rel_len,
		"stats/assess_assembly_all_meanQ.tsv",
		"stats/assess_homopolymers_all_rel_len.tsv",
		"stats/assess_homopolymers_all_correct_len.tsv",
		list_assess_homopolymers_correct_len_pdf,
		list_busco_out

rule assess_assembly:
	input:
		assembly = "assemblies/{id}.fa",
		ref = "references/{ref}.fa"
	output:
		summ = "{id}/assess_assembly/{id}_{ref}_summ.txt",
		meanQ = "{id}/assess_assembly/{id}_{ref}_meanQ.tsv"
	log: "{id}/assess_assembly/{id}_{ref}_log.txt"
	shell:
		"""
		assess_assembly -r {input.ref} -i {input.assembly} -p {wildcards.id}/assess_assembly/{wildcards.id}_{wildcards.ref} >{log} 2>&1
		meanQ=$(grep -A2 "{wildcards.ref} Q" {output.summ} | tail -n1 | awk '{{print $2}}')
		echo "{wildcards.id}\t{wildcards.ref}\t$meanQ" > {output.meanQ}
		"""

rule assess_homopolymers_minimap:
	threads: 5
	input:
		assembly = "assemblies/{id}.fa",
		ref = "references/{ref}.fa"
	output: "{id}/assess_homopolymers/{id}_{ref}.bam"
	log: "{id}/assess_homopolymers/{id}_{ref}.minimap2.log"
	shell:
		"""
		minimap2 -x asm5 -t {threads} --MD -a {input.ref} {input.assembly} 2>{log} | samtools sort -o {output} - >>{log} 2>&1
		"""

rule assess_homopolymers:
	input: "{id}/assess_homopolymers/{id}_{ref}.bam"
	output:
		count = "{id}/assess_homopolymers/{id}_{ref}_count/hp_counts.pkl",
		rel_len = "{id}/assess_homopolymers/{id}_{ref}_analyse/hp_rel_len_counts.txt",
		rel_len2 = "{id}/assess_homopolymers/{id}_{ref}_rel_len.tsv",
		correct_len = "{id}/assess_homopolymers/{id}_{ref}_analyse/hp_correct_vs_len.txt",
		correct_len2 = "{id}/assess_homopolymers/{id}_{ref}_correct_len.tsv"
	params:
		count_dir = directory("{id}/assess_homopolymers/{id}_{ref}_count"),
		analyse_dir = directory("{id}/assess_homopolymers/{id}_{ref}_analyse")
	log: "{id}/assess_homopolymers/{id}_{ref}_log.txt"
	shell:
		"""
		rm -rf {params.count_dir} {params.analyse_dir}
		assess_homopolymers count -o {params.count_dir} {input} >{log} 2>&1
		assess_homopolymers analyse -o {params.analyse_dir} {output.count} >>{log} 2>&1

		grep -v count {output.rel_len} | sed 's/^/{wildcards.id}\t{wildcards.ref}\t/' > {output.rel_len2}
		grep -v rlen {output.correct_len} | sed 's/^/{wildcards.id}\t{wildcards.ref}\t/' > {output.correct_len2}
		"""

rule gather_stats:
	input:
		aa_meanQ = list_assess_assembly_meanQ,
		hp_rel_len = list_assess_homopolymers_rel_len,
		hp_correct_len = list_assess_homopolymers_correct_len,
	output:
		all_aa_meanQ = "stats/assess_assembly_all_meanQ.tsv",
		all_hp_rel_len = "stats/assess_homopolymers_all_rel_len.tsv",
		all_hp_correct_len = "stats/assess_homopolymers_all_correct_len.tsv",
	shell:
		"""
		#filter inf values, which happen when comparing identical assemblies
		cat {input.aa_meanQ} | grep -v -w 'inf' > {output.all_aa_meanQ}

		cat {input.hp_rel_len}  > {output.all_hp_rel_len}
		cat {input.hp_correct_len}  > {output.all_hp_correct_len}
		"""

rule plot:
	input:
		all_aa_meanQ = "stats/assess_assembly_all_meanQ.tsv",
		all_hp_correct_len = "stats/assess_homopolymers_all_correct_len.tsv",
	output:
		all_correct_len_pdfs = list_assess_homopolymers_correct_len_pdf
	script: "scripts/plot.R"

rule busco:
	threads: 5
	input:
		assembly = "assemblies/{id}.fa"
	output:
		"{id}/busco/short_summary.specific.{busco_lineage}_odb10.busco.txt"
	log: "{id}/busco/busco_{busco_lineage}.log"
	shell:
		"""
		cd {wildcards.id} && busco -c {threads} -f -m genome -l {busco_lineage} -o busco -i ../{input} >../{log} 2>&1
		"""
