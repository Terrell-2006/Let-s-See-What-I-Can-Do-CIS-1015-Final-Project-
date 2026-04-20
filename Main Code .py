# Validate DNA 
def validate_dna(seq):
    return all(base in "ATCG" for base in seq)

# Transcription 
def transcribe(seq):
    return seq.replace("T", "U")

# --- Translation ---
def translate(rna):
    protein = ""
    for i in range(0, len(rna)-2, 3):
        codon = rna[i:i+3]
        amino = codon_table.get(codon, "?")
        if amino == "STOP":
            break
        protein += amino
    return protein

# --- GC Content ---
def gc_content(seq):
    return round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2)
