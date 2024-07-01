def read_gene_ids(file_path):
    """Read gene IDs from a text file."""
    gene_ids = set()
    with open(file_path, 'r') as file:
        for line in file:
            gene_id = line.strip()
            gene_ids.add(gene_id)
    return gene_ids

def extract_sequences(fasta_file, genes_of_interest_file, output_file):
    """Extract sequences for genes of interest from a FASTA file."""
    gene_ids_of_interest = read_gene_ids(genes_of_interest_file)
    
    current_gene_id = None
    include_current_gene = False

    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as output:
        for line in fasta:
            if line.startswith('>'):
                current_gene_id = line.strip()[1:]  # Remove '>'
                include_current_gene = current_gene_id in gene_ids_of_interest
            if include_current_gene:
                output.write(line)

# define input files, 
# all_genes_fasta.fa is the fasta file with all gene sequences extracted from it ( i saved this specific file to misc folder)
# genes_of_interest_file.txt includes a text file with a list of the gene names 
# specify the output file with output.fasta
fasta_file = 'all_genes_fasta.fa'
genes_of_interest_file = 'genes_of_interest.txt'
output_file = 'genes_of_interest.fasta'

extract_sequences(fasta_file, genes_of_interest_file, output_file)
