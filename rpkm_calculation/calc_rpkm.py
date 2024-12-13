import sys

def calc_rpkm(reads: str, gene_lengths: dict, total_reads: int):
    """
    :reads: Sorted and filtered .annot file
    :gene_lengths: A dict containing gene ids and lengths
    """

    per_million_scaling_factor = total_reads / 1000000
    with open(reads) as read_file:

        read_entries = read_file.readlines()
        read_counter_per_gene = 0
        last_gene_name = "-1"

        for entry in read_entries:
            main_comp = entry.split("\t")
            gene_name = main_comp[1].split(",")[0]
            if gene_name == last_gene_name:
                read_counter_per_gene += 1
            elif (last_gene_name == "-1"):  # first ieration
                last_gene_name = gene_name
                read_counter_per_gene += 1
            else:
                RPKM = read_counter_per_gene / ( (gene_lengths[last_gene_name][0] / 1000) * per_million_scaling_factor)
                print(f"{last_gene_name}\t{RPKM}\t{gene_lengths[last_gene_name][1].strip("\n")}")
                # print(f"{last_gene_name}: RPKM = {read_counter_per_gene} / (({gene_lengths[last_gene_name][0]} / 1000) * {per_million_scaling_factor}) = {RPKM}")
                
                read_counter_per_gene = 1  # reset read count
                last_gene_name = gene_name
        
        if last_gene_name in gene_lengths:
            RPKM = read_counter_per_gene / ( (gene_lengths[last_gene_name][0] / 1000) * per_million_scaling_factor)
            print(f"{last_gene_name}\t{RPKM}\t{gene_lengths[last_gene_name][1].strip("\n")}")
            # print(f"{last_gene_name}: RPKM = {read_counter_per_gene} / ( {gene_lengths[last_gene_name][0]} / 1000 * {per_million_scaling_factor})")


def read_gene_lengths(gene_path):
    gene_lengths = dict()
    with open(gene_path) as file:
        lines = file.readlines()
        for line in lines:
            if(line.startswith("gene")):
                continue
            comp = line.split("\t")
            gene_lengths[comp[0]] = [int(comp[1]), comp[2]]

    return gene_lengths
 
if __name__ == "__main__":
    path_to_reads = sys.argv[1]
    path_to_gene_lengths = sys.argv[2]
    total_reads = int(sys.argv[3])
    gene_lengths = read_gene_lengths(path_to_gene_lengths)
    calc_rpkm(path_to_reads, gene_lengths, total_reads)
