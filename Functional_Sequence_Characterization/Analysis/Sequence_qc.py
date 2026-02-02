from Bio import SeqIO

record = SeqIO.read("Input_sequence.fasta", "fasta")

sequence = record.seq
length = len(sequence)

print("-" * 120)
print("Protein ID:", record.id)
print("Description:", record.description)
print("Protein Length:", length,"aa")
print("Sequence:", sequence)
print("-" * 120)

aa_weights = {
    'A': 89.09,  'R': 174.20, 'N': 132.12, 'D': 133.10,
    'C': 121.15, 'E': 147.13, 'Q': 146.15, 'G': 75.07,
    'H': 155.16, 'I': 131.17, 'L': 131.17, 'K': 146.19,
    'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09,
    'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}

molecular_weight = 0

for aa in sequence:
    if aa in aa_weights:
        molecular_weight += aa_weights[aa]

water_weight = 18.015
molecular_weight = molecular_weight - (water_weight * (length - 1))

print("Estimated Molecular Weight: %.2f Da" % molecular_weight)
print("-" * 120)

amino_acids = "ACDEFGHIKLMNPQRSTVWY"
print("Amino Acid Composition (%):")
for aa in amino_acids:
    percentage = sequence.count(aa) / length * 100
    print("%s: %.2f%%" % (aa, percentage))

print("-" * 120)

basic = "KRH"
acidic = "DE"
aromatic = "FWY"
polar = "STNQC"
non_polar = "AVILMG"

basic_count = 0
acidic_count = 0
aromatic_count = 0
polar_count = 0
non_polar_count = 0

for aa in sequence:
    if aa in basic:
        basic_count += 1
    if aa in acidic:
        acidic_count += 1
    if aa in aromatic:
        aromatic_count += 1
    if aa in polar:
        polar_count += 1
    if aa in non_polar:
        non_polar_count += 1

print("Amino Acid Group Composition (%):")
print("Basic: %.2f%%" % ((basic_count / length) * 100))
print("Acidic: %.2f%%" % ((acidic_count / length) * 100))
print("Aromatic: %.2f%%" % ((aromatic_count / length) * 100))
print("Polar: %.2f%%" % ((polar_count / length) * 100))
print("Non-polar: %.2f%%" % ((non_polar_count / length) * 100))
print("-" * 120)

if len(sequence) < 100:
    print("Note: Sequence too short for further analysis.")
else:
    print("Note: Sequence accepted for downstream analysis.")
print("-" * 120)