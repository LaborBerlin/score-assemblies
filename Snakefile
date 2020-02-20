shell.executable("/bin/bash")

assemblies, = glob_wildcards("assemblies/{id}.fa")
references, = glob_wildcards("references/{id}.fa")


rule all:
	input:
		expand("{id}/assess_assembly/{id}_{ref}_stats.txt", id=assemblies, ref=references),
		expand("{id}/assess_homopolymers/{id}_{ref}.bam", id=assemblies, ref=references),
		expand("{id}/assess_homopolymers/{id}_{ref}_analyse/hp_correct_vs_len.txt", id=assemblies, ref=references)

rule assess_assembly:
	input:
		assembly = "assemblies/{id}.fa",
		ref = "references/{ref}.fa"
	output: "{id}/assess_assembly/{id}_{ref}_stats.txt"
	log: "{id}/assess_assembly/{id}_{ref}_log.txt"
	shell:
		"""
		assess_assembly -r {input.ref} -i {input.assembly} -p {wildcards.id}/assess_assembly/{wildcards.id}_{wildcards.ref} >{log} 2>&1
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
		analyse = "{id}/assess_homopolymers/{id}_{ref}_analyse/hp_correct_vs_len.txt"
	params:
		count_dir = directory("{id}/assess_homopolymers/{id}_{ref}_count"),
		analyse_dir = directory("{id}/assess_homopolymers/{id}_{ref}_analyse")
	log: "{id}/assess_homopolymers/{id}_{ref}_log.txt"
	shell:
		"""
		rm -rf {params.count_dir} {params.analyse_dir}
		assess_homopolymers count -o {params.count_dir} {input} >{log} 2>&1
		assess_homopolymers analyse -o {params.analyse_dir} {output.count} >>{log} 2>&1
		"""

