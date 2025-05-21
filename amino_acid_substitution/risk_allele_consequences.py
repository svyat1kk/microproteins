from pyfaidx import Fasta
import pandas as pd
from Bio.Seq import Seq

genome = Fasta("Human_Genome_Full_GRCh38.fa")
remainings = pd.read_csv("./remainings_to_analyse_14.txt", sep = '\t')

remainings["CHR_ID"] = remainings["CHR_ID"].apply(lambda x: f"chr{x}")

aa_groups = {
    "Aliphatic": set(["G", "A", "V", "L", "I"]),
    "Hydroxyl/Sulfur": set(["S", "C", "U", "T", "M"]),
    "Cyclic": set(["P"]),
    "Aromatic": set(["F", "Y", "W"]),
    "Basic": set(["H", "K", "R"]),
    "Acidic/Amide": set(["D", "E", "N", "Q"])
}

def get_aa_group(aa):
    for group, residues in aa_groups.items():
        if aa in residues:
            return group
    return "unknown group"


def snp_risk_allele(seq, chr_id, snp_pos,risk_allele):
    snp_rel_pos = int(float(snp_pos)) - int(seq['start'])
    seq_str = str(seq['sequence'])
    mutated = seq_str[:snp_rel_pos] + risk_allele + seq_str[snp_rel_pos+1:]
    return mutated

def convert_to_protein(dna_seq):
        try:
            if dna_seq[:3] == "ATG":
                wild_protein = Seq(dna_seq).translate()
            else:
                reverse_complement = Seq(dna_seq).reverse_complement()
                wild_protein = reverse_complement.translate()
        except:
            print("Could not translate the DNA sequence")
            wild_protein = "Could not be translated"

        return wild_protein

def point_out_mutation(wild, mutated):
    for i, aa in enumerate(wild):
        if wild[i] != mutated[i]:
            wild_group = get_aa_group(wild[i])
            mutated_group = get_aa_group(mutated[i])
            message = f"{i+1}:{wild[i]}->{mutated[i]} | Group: {wild_group} → {mutated_group}"
            break
        else:
            pass
    
    return message

results = []

for _, row in remainings.iterrows():
    chrom = str(row['CHR_ID'])
    start = int(row['starts']) - 1
    end = int(row['ends'])

    snp_pos_str = row["CHROM"]
    snp_chr, snp_pos = snp_pos_str.split("_")
    risk_allele = row["STRONGEST SNP-RISK ALLELE"].split("-")[1]

    seq = genome[chrom][start:end].seq  
    wild_seq = str(seq)

    rel_pos = int(float(snp_pos)) - start
    if 0 <= rel_pos < len(wild_seq):
        mutated_seq = wild_seq[:rel_pos] + risk_allele + wild_seq[rel_pos + 1:]
    else:
        mutated_seq = None

    wild_protein = convert_to_protein(wild_seq)
    mutated_protein = convert_to_protein(mutated_seq)

    if wild_protein != mutated_protein:
        mutation = point_out_mutation(wild_protein, mutated_protein)
    else:
        mutation = "did not occur"

    results.append({
        "chr_id": chrom,
        "start": start,
        "end": end,
        "snp_id": row["STRONGEST SNP-RISK ALLELE"],
        "snp_position": snp_pos_str,
        "orf_name": row["orf_name"],
        "aminoacid_seq_mutation": mutation,
        "wild_protein": wild_protein,
        "mutated_protein": mutated_protein,
        "wild_ORF_seq": wild_seq,
        "mutated_ORF_seq": mutated_seq
    })

    
output_df = pd.DataFrame(results)
for i in range(remainings.shape[0]):
    if remainings.orf_sequence[i] != output_df.wild_protein[i]:
        print(f"The protein sequence obtained from reference human genome does not match with orf_sequence in row №{i+1} from the dataset with microproteins")
output_df.to_csv("./mutated_sequences_14.txt", sep='\t', index=False)