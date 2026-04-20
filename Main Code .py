# DNA Analyzer & Mutation Classifier

# Full RNA Codon Table 
codon_table = {
    'UUU':'F','UUC':'F','UUA':'L','UUG':'L',
    'CUU':'L','CUC':'L','CUA':'L','CUG':'L',
    'AUU':'I','AUC':'I','AUA':'I','AUG':'M',
    'GUU':'V','GUC':'V','GUA':'V','GUG':'V',
    'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
    'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
    'UAU':'Y','UAC':'Y','UAA':'STOP','UAG':'STOP',
    'CAU':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAU':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAU':'D','GAC':'D','GAA':'E','GAG':'E',
    'UGU':'C','UGC':'C','UGA':'STOP','UGG':'W',
    'CGU':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGU':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGU':'G','GGC':'G','GGA':'G','GGG':'G'
}

# Validate DNA 
def validate_dna(seq):
    return all(base in "ATCG" for base in seq)

# Transcription 
def transcribe(seq):
    return seq.replace("T", "U")

# Translation
def translate(rna):
    protein = ""
    for i in range(0, len(rna)-2, 3):
        codon = rna[i:i+3]
        amino = codon_table.get(codon, "?")
        if amino == "STOP":
            break
        protein += amino
    return protein

# GC Content
def gc_content(seq):
    return round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2)

# Mutation Classification (Fix)
def classify_mutation(seq1, seq2):
    mutations = []   
    
    for i in range(0, min(len(seq1), len(seq2)), 3):  
        codon1 = transcribe(seq1[i:i+3])   
        codon2 = transcribe(seq2[i:i+3])
        
        if codon1 != codon2 and len(codon1) == 3:   
            aa1 = codon_table.get[codon1, "?"]   
            aa2 = codon_table.get(codon2, "?")  
            
            if aa1 = aa2:   
                mtype == "Silent"   
            elif aa2 == "STOP":   
                mtype = "Nonsense"  
            else:
                mtype = "Missense"   
            
            mutations.append((i, codon1, codon2, aa1, aa2, mtype))  
    
    return mutations

# FASTA file reader (Fix)
def readFasta(filePath):
    seq = ""   
    file = open(filePath, "r")   
    for line in file:
        if line[0] != ">":   
            seq += line.strip   
        else:
            continue     
        if line == "\n":   
            pass
    
    file.close   
    
return seq.upper    
