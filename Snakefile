#this workflow reconstructs the mutations in the RSV A and RSV B duplicated regions of the G gene

A_OR_B = ["a"]

rule all:
    input:
        expand("results/{a_or_b}_indel_with_effect.tsv", a_or_b=A_OR_B)

rule branch_from_root:
    message:
        """Adding mutations to root sequence to reconstruct all branches"""
    input:
        reference = "data/{a_or_b}reference.gbk",
        root = "data/{a_or_b}/root_sequence.json",
        sequences = "data/{a_or_b}/alignment.fasta",
        tree = "data/{a_or_b}/tree.nwk",
        treejson = "data/rsv_{a_or_b}_genome.json"
    params:
        gene = "G"
    output:
        reconstructed_seq = "results/{a_or_b}/reconstructed_gene.fasta"
    shell:
        """
        python3 scripts/reconstruct_sequences.py \
        --root {input.root} \
        --json {input.treejson} \
        --alignment {input.sequences} \
        --tree {input.tree} \
        --output {output.reconstructed_seq} \
        --gene {params.gene} \
        --reference {input.reference}
        """


rule insertions_deletions:
    input:
        reference ="data/{a_or_b}reference.gbk",
        tree_ = "data/{a_or_b}/tree.nwk",
        insertions = "data/{a_or_b}/insertions.csv",
        jsonfile = "data/rsv_{a_or_b}_genome.json"
    output:
        tsv = "results/{a_or_b}_indel.tsv"
    shell:
        """
        python3 scripts/insertions_and_deletions.py \
        --json {input.jsonfile} \
        --tree {input.tree_} \
        --reference {input.reference} \
        --output {output.tsv} \
        --insertions {input.insertions}
        """

rule compare:
    input:
        G_aligned = "data/{a_or_b}/alignedG.fasta",
        indels = rules.insertions_deletions.output
    output:
        "results/{a_or_b}_indel_with_effect.tsv"
    shell:
        """
        python3 scripts/compare.py \
        --alignedgene {input.G_aligned} \
        --indels {input.indels} \
        --output {output}
        """

    