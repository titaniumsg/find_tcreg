#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: May 27, 2024
#   Email: sagarg@sas.upenn.edu

# Import required libraries
import csv
from collections import defaultdict
import argparse
from pymol import cmd, stored
from pathlib import Path

from Bio.Align import substitution_matrices
blosum62_matrix = substitution_matrices.load("BLOSUM62")
pos_charge = ['R', 'H', 'K']
neg_charge = ['D', 'E', 'Q', 'N', 'S', 'T']

# Import custom libraries
from fasta import FASTA

CUTOFF_DISTANCE = 5  # Angstrom cutoff to define a "contact"

def parse_args():

    parser = argparse.ArgumentParser(description='Finds T-CREGs given a list of HLAs that bind to a peptide')
    parser.add_argument("-hla_file", help="File containing HLAs that bind to the peptide; needs to follow the format A*02:01", required=True)
    parser.add_argument("-output_dir", help="Output directory", required=True)

    parser.add_argument("-tcreg_file", help="Starting T-CREG file. Default is the 21 T-CREG, >5 percent cutoff version. Can be switched to tcregs_greedy_no_cutoff.txt", default="tcregs_greedy_5_percent_cutoff.txt")

    # If you have a structure, you can use it to get more specific TCR contact residues
    parser.add_argument("-structure", help="File containing pHLA-receptor structure", type=str, default="")
    parser.add_argument("-hla_chain", help="Chain name of HLA in the pHLA-receptor structure", type=str, default="A")
    parser.add_argument("-receptor_chain", nargs="+", help="Chain name(s) of the receptor in the pHLA-receptor structure. If you have multiple chain names, give them as -receptor_chain D E F", default=["D", "E"])
    return parser.parse_args()

def get_mutation_list(seq1, seq2):

    seq1_mut, seq2_mut = [], []
    for seq1_aa, seq2_aa in zip(seq1, seq2):
        if seq1_aa != seq2_aa:
            seq1_mut.append(seq1_aa)
            seq2_mut.append(seq2_aa)
    
    return seq1_mut, seq2_mut

def find_tcregs(seq_tcr_contact, new_alleles):
        
    seq_max = {}
    for allele1 in new_alleles:
        seq1 = seq_tcr_contact[allele1]
        seq_max[allele1] = [allele1]
        for allele2 in new_alleles:
            if allele1 == allele2: continue
            
            seq2 = seq_tcr_contact[allele2]
            
            # If sequences are identical, group them together
            if seq1 == seq2: 
                seq_max[allele1].append(allele2)
                if allele1 not in seq_max[allele1]:
                    seq_max[allele1].append(allele1)
                continue
            
            # Check if mutations between sequences are significant
            seq1_mut, seq2_mut = get_mutation_list(seq1, seq2)
            significant_mut = False
            
            for aa1, aa2 in zip(seq1_mut, seq2_mut):
                try:
                    blosum_score = blosum62_matrix[(aa1, aa2)]
                except:
                    blosum_score = blosum62_matrix[(aa2, aa1)]

                # Check for charge switches (see Table S1)
                charge_switch = (
                    (aa1 in pos_charge and aa2 in neg_charge) or 
                    (aa1 in neg_charge and aa2 in pos_charge)
                )

                if blosum_score < 0 or charge_switch:
                    significant_mut = True
                    break
        
            if significant_mut:
                continue

            else:
                seq_max[allele1].append(allele2)
                if allele1 not in seq_max[allele1]:
                    seq_max[allele1].append(allele1)
    
    return seq_max

def write_tcreg_file(filename, tcreg, usa_allelic_freq, tcr_contact_resis, seq_tcr_contact):

    def get_mutations_dict(sequences, allele_list):

        mutations = defaultdict(set)
        for allele1 in allele_list:
            seq1 = sequences[allele1]
            for allele2 in allele_list:
                seq2 = sequences[allele2]
                if seq1 == seq2 or allele1 == allele2: continue

                for seq1_aa, seq2_aa, hla_pos in zip(seq1, seq2, tcr_contact_resis):
                    if seq1_aa != seq2_aa:
                        allele_combo = sorted([allele1, allele2])
                        wt_aa = sequences[allele_combo[0]][tcr_contact_resis.index(hla_pos)]
                        mut_aa = sequences[allele_combo[1]][tcr_contact_resis.index(hla_pos)]
                        mutations["_".join(allele_combo)].add(f"{wt_aa}{hla_pos}{mut_aa}")
                        
        return mutations

    final_file = open(filename, "w")
    
    for common_allele, allele_list in tcreg.items():
        total_allelic_coverage = 0
        for allele in allele_list:
            total_allelic_coverage += usa_allelic_freq[allele]

        # Get unique mutations within this T-CREG
        mutations = []
        all_mutations_dict = get_mutations_dict(seq_tcr_contact, allele_list)
        for allele_pair, mut_set in all_mutations_dict.items():
            mut_list = list(mut_set)
            for mut in mut_list:
                pos = mut[1:-1]
                aa1 = mut[0]
                aa2 = mut[-1]
                if mut in mutations or f"{aa2}{pos}{aa1}" in mutations: 
                    continue
                else:
                    mutations.append(mut)
        
        print(f"{common_allele} | {total_allelic_coverage*100:.2f}% | {', '.join(mutations)}", file=final_file)
        for allele in allele_list:
            print(f" - {allele}", file=final_file)

def get_tcreg(hla_file, output_dir, tcreg_file, structure, hla_chain, receptor_chain):

    '''
    Takes as input a peptide sequence and list of HLA alleles that are assumed to be binders

    Finds overlap between strong binders and T-CREG alleles
    Creates sequence logos of non-conserved TCR contact residues for strong binders
    '''

    # Consider HLA (-A, -B, -C) alleles that have an allelic frequency of ≥0.05% in at least one of the 5 major American sub-populations (n = 215)
    file = FASTA("imgt_allelic_freq_0.05_NMDP.fasta")  # IPD-IMGT/HLA Release 3.48.0
    file.read()
    hla_seq = file.get_sequences()

    # Load population and allelic frequency values for each allele from the NMDP/BeTheMatch database
    usa_allelic_freq = defaultdict(float)
    with open("NMDP_per_allele.csv", "r") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for line in reader:
            usa_allelic_freq[line[0]] = float(line[-2])

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if structure == "":
        # Use pre-computed TCR contact residues from TCR3d database
        tcr_contact_resis = []
        with open(f"tcr_contact_resi_{CUTOFF_DISTANCE}.csv", "r") as csvfile:
            reader = csv.reader(csvfile)
            next(reader)
            for line in reader:
                # Only include residues that are contacts in >10% of structures
                if float(line[-1]) > 10:
                    tcr_contact_resis.append(int(line[1]))
        tcr_contact_resis = sorted(tcr_contact_resis)
        
        print(f"Evaluating T-CREGs using {len(tcr_contact_resis)} TCR contact residues ({CUTOFF_DISTANCE} A cutoff) based on structures in TCR3d")

        # Load existing T-CREG definitions from file
        tcreg = defaultdict(list)
        with open(tcreg_file, "r") as txtfile:
            lines = txtfile.readlines()
            lines = [l.strip() for l in lines]
            current_tcreg = ""   
            for line in lines:
                if "-" in line:
                    tcreg[current_tcreg].append(line.split(" ")[-1])
                elif "|" in line:
                    current_tcreg = line.split(" ")[0]
                    tcreg[current_tcreg] = []

    else:
        # Use provided structure to determine TCR contact residues
        if not Path(structure).is_file():
            raise FileNotFoundError(f"{structure} could not be found! Was the correct path to the structure given?")
        
        receptor_chain_str = " or ".join([f"chain {chain}" for chain in receptor_chain])
        
        stored.tcr_contact_resis = []
        cmd.load(structure, "struc")
        cmd.remove("hydro and struc") # remove hydrogen atoms
        cmd.select("interface", f"(br. chain {hla_chain} within {CUTOFF_DISTANCE} of ({receptor_chain_str})) and (chain {hla_chain} and (sc. or name ca))")
        cmd.iterate("interface and name ca", "stored.tcr_contact_resis.append(resi)")

        if len(stored.tcr_contact_resis) == 0:
            raise ValueError(f"Could not find any TCR contact residues in {structure}. Are the chain names correct?")

        tcr_contact_resis = sorted([int(pos) for pos in stored.tcr_contact_resis])
        print(f"Generating T-CREGs using {len(tcr_contact_resis)} TCR contact residues ({CUTOFF_DISTANCE} A cutoff) based off of {structure}")
        cmd.delete("all")
        cmd.quit()
        
        # Extract TCR contact residue sequences for each allele
        seq_tcr_contact = defaultdict(str)
        for allele, seq in hla_seq.items():
            str_tcr_contact = ''
            for index, aa in enumerate(list(seq)):
                pos = index+1
                if pos in tcr_contact_resis: # only consider tcr contact positions
                    str_tcr_contact += aa
            seq_tcr_contact[allele] = str_tcr_contact

        seq_max = {}
        seq_max = find_tcregs(seq_tcr_contact, list(seq_tcr_contact.keys()))
        new_alleles = list(seq_tcr_contact.keys())

        # Greedy algorithm to determine T-CREGss
        temp_tcreg = {}
        while True:

            seq_max = find_tcregs(seq_tcr_contact, new_alleles)

            seq_length = defaultdict(float)
            for allele, a_list in seq_max.items():
                seq_length[allele] = len(a_list)
            
            seq_length = {k: v for k, v in sorted(seq_length.items(), key=lambda item: item[1], reverse=True)}

            # Temporary representative - this will be changed to the most common allele in the T-CREG later
            best_rep = list(seq_length.keys())[0]
            temp_tcreg[best_rep] = seq_max[best_rep]

            # Remove the alleles that we have covered with the identified T-CREG
            for allele in seq_max[best_rep]:
                new_alleles.remove(allele)
            for allele in new_alleles:
                del seq_max[allele]
            
            # print(f"Taking out {best_rep} and its {len(seq_max)} neighbors leaving {len(new_alleles)} alleles")
            if len(new_alleles) == 0: break

        # Convert temporary T-CREGs to final format with most common allele as representative
        temp_tcreg = dict(sorted(temp_tcreg.items(), key=lambda x:x[0]))
        tcreg = defaultdict(list)
        for temp_allele, allele_list in temp_tcreg.items():
            total_allelic_coverage = 0
            allele_freq = {}
            
            for allele in allele_list:
                
                total_allelic_coverage += usa_allelic_freq[allele]
                allele_freq[allele] = usa_allelic_freq[allele]

                allele_freq = dict(sorted(allele_freq.items(), key=lambda x:x[1], reverse=True))
            
                # if total_allelic_coverage > 0.05:
            most_common_allele = list(allele_freq.keys())[0]
            
            tcreg[most_common_allele] = allele_list
        
        filename = f"{output_dir}/{Path(structure).stem}_tcreg.txt"
        write_tcreg_file(filename, tcreg, usa_allelic_freq, tcr_contact_resis, seq_tcr_contact)

        print(f"A new T-CREG list was generated in {filename}")
        

    with open(hla_file) as txtfile:
        binders = txtfile.readlines()
        binders = [allele.strip().replace("HLA", "").replace("-", "") for allele in binders]
    
    # Extract TCR contact sequences for all HLA alleles
    seq_tcr_contact = defaultdict(str)
    for allele, seq in hla_seq.items():
        str_tcr_contact = ''
        for index, aa in enumerate(list(seq)):
            pos = index+1
            if pos in tcr_contact_resis: # only consider tcr contact positions
                str_tcr_contact += aa
        seq_tcr_contact[allele] = str_tcr_contact

    tcreg_binders = defaultdict(list)
    considered_alleles = []
    for allele in binders:
        for creg, a_list in tcreg.items():
            if allele in a_list:
                tcreg_binders[creg].append(allele)
                considered_alleles.append(allele)
    
    if len(tcreg_binders) == 0:
        raise ValueError(f"No T-CREGs found that overlap with binders. Is your HLA list ({hla_file}) correct?")
    
    with open(f"{output_dir}/tcreg_final.txt", "w") as txtfile:
        for creg, allele_list in tcreg_binders.items():
            
            total_allelic_freq = 0
            for allele in allele_list:
                total_allelic_freq += usa_allelic_freq[allele]

            txtfile.write(f"{creg} T-CREG | {total_allelic_freq*100:.2f}%\n")
            for allele in allele_list:
                txtfile.write(f" - {allele} | {usa_allelic_freq[allele]*100:.2f}%\n")
    
if __name__ == "__main__":

    args = parse_args()
    get_tcreg(args.hla_file, args.output_dir, args.tcreg_file, args.structure, args.hla_chain, args.receptor_chain)
