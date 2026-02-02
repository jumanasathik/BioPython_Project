from Bio.Blast import NCBIWWW
from Bio import SeqIO

record = SeqIO.read("Input_sequence.fasta", "fasta")
print("-" * 120)
print("Running BLASTp... please wait (this may take time)")

result_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence=record.seq
)

with open("blast_result.xml", "w") as bp:
    bp.write(result_handle.read())

print("BLASTp search completed successfully!")
print("Results saved as blast_result.xml")
print("-" * 120)

from Bio.Blast import NCBIXML

with open("blast_result.xml") as b:
    blast_record = NCBIXML.read(b)

high_identity_hits = []
for hit in blast_record.alignments:
    for hsp in hit.hsps:
        identity = (hsp.identities / hsp.align_length) * 100
        if identity >= 90:
            high_identity_hits.append({
                "hit_id": hit.hit_id,
                "hit_def": hit.hit_def,
                "length": hit.length,
                "identity": identity,
                "evalue": hsp.expect,
                "query_start": hsp.query_start,
                "query_end": hsp.query_end,
                "sbjct_start": hsp.sbjct_start,
                "sbjct_end": hsp.sbjct_end,
                "query_aln": hsp.query,
                "subject_aln": hsp.sbjct,
                "match_line": hsp.match
            })
            break  

print("Total hits with >=90% dentity:", len(high_identity_hits))
print("-" * 120)

print("\n----------------------------- Closest Homologs -----------------------------")
for hit in high_identity_hits:
    print("\nID:", hit["hit_id"])
    print("Definition:", hit["hit_def"])
    print("Length:", hit["length"])
    print("Percent Identity: %.2f%%" % hit["identity"])
    print("E-value:", hit["evalue"])
print("-" * 120)

print("\n----------------------------- Alignments & Conserved Regions -----------------------------")
for hit in high_identity_hits:
    print("\nID:", hit["hit_id"])
    print("Query Alignment Positions: %d to %d" % (hit["query_start"], hit["query_end"]))
    print("Subject Alignment Positions: %d to %d" % (hit["sbjct_start"], hit["sbjct_end"]))
    print("Conserved Region:")
    print(hit["match_line"])
print("-" * 120)

print("\n----------------------------- Evolutionary Hints -----------------------------")
for hit in high_identity_hits:
        identity = hit["identity"]
        evalue = hit["evalue"]

        print("\nID: " + hit["hit_id"] + "\n")

        if identity == 100:
            print("Sequences are identical. Likely very close evolutionary relationship.\n")
        elif identity >= 99:
            print("Almost identical, likely same or very recent homolog.\n")
        elif identity >= 90:
            print("Very close homolog, strong evolutionary relationship.\n")
        elif identity >= 70:
            print("Moderate similarity, possible distant homolog.\n")
        else:
            print("Low similarity, likely distant relationship.\n")

        if evalue < 1e-50:
            print("Extremely significant match (very reliable homology).\n")
        elif evalue < 1e-10:
            print("Highly significant match.\n")
        else:
           print("Moderate significance.\n")

print("-" * 120 + "\n")

print("\n----------------------------- Functional Annotation -----------------------------")

if len(high_identity_hits) > 0:
    top_hit = high_identity_hits[0]
    description = top_hit["hit_def"]

    print("\nTop Closest Homolog:")
    print("Definition:", description)
    print("Percent Identity: %.2f%%" % top_hit["identity"])
    print("E-value:", top_hit["evalue"])

    print("\nPredicted Function:")
    print("The BLAST hit annotation suggests that this protein is related to:")
    print(description)
    print("\nSince the query sequence shows high identity and a very low E-value,")
    print("it is likely that the query protein performs a similar biological function.")
else:
    print("No high-identity homologs found. Functional prediction is uncertain.")

print("-" * 120)

