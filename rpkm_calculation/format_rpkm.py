import sys

def split_genes(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip("\n")
            if "|" in line:
                main_comp = line.split("\t")
                id = main_comp[0]
                pcr = main_comp[-1]
                genes = main_comp[1].split("|")
                for g in genes:
                    print(f"{id}\t{g}\t{pcr}")
            else:
                print(line)

if __name__ == "__main__":
    path = sys.argv[1]
    split_genes(path)
