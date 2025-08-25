import random

def generate_artificial_bed(
    genome_length=1_000_000,   # length of chromosome
    num_genes=50,               # total number of genes
    min_gene_length=1000,       # minimum gene length
    max_gene_length=100000,     # maximum gene length
    chromosomes=["chr1"],       # list of chromosomes
    output_file="artificial_genes.bed"
):
    # generate some artificial gene names
    gene_prefixes = ["GEN", "ART", "LOC", "FAKE", "TEST"]
    
    with open(output_file, "w") as f:
        for i in range(num_genes):
            chrom = random.choice(chromosomes)
            gene_length = random.randint(min_gene_length, max_gene_length)
            start = random.randint(0, genome_length - 1)
            end = start + gene_length
            gene_name = f"{random.choice(gene_prefixes)}{i+1}"
            f.write(f"{chrom}\t{start}\t{end}\t{gene_name}\n")

    print(f"Generated {num_genes} artificial genes in {output_file}")

# Example usage:
generate_artificial_bed(
    genome_length=2_500_000,
    num_genes=100,
    chromosomes=["chr1"], #chromosomes=["chr1","chr2","chr3"],
    min_gene_length=5000,
    max_gene_length=50000,
    output_file="artificial_overlapping_genes.bed"
)
