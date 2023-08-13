#this workflow reconstructs the mutations in the RSV A and RSV B duplicated regions of the G gene

A_OR_B = ["a", "b"]

rule all:
    input:
        expand("results/{a_or_b}_indel.tsv", a_or_b=A_OR_B)

rule insertions_deletions:
    input:
        reference ="data/{a_or_b}reference.gbk",
        tree_ = "data/{a_or_b}/tree.nwk",
        insertions = "data/{a_or_b}/insertions.csv",
        jsonfile = "data/rsv_{a_or_b}_genome.json",
        G_aligned = "data/{a_or_b}/alignedG.fasta",
    output:
        tsv = "results/{a_or_b}_indel.tsv"
    shell:
        """
        python3 scripts/insertions_and_deletions.py \
        --json {input.jsonfile} \
        --tree {input.tree_} \
        --reference {input.reference} \
        --output {output.tsv} \
        --insertions {input.insertions} \
        --alignedgene {input.G_aligned}
        """

#rule compare:
#    input:
#        G_aligned = "data/{a_or_b}/alignedG.fasta",
#        indels = rules.insertions_deletions.output
#    output:
#        "results/{a_or_b}_indel_with_effect.tsv"
#    shell:
#        """
#        python3 scripts/compare.py \
#        --alignedgene {input.G_aligned} \
#        --indels {input.indels} \
#        --output {output}
#        """

    