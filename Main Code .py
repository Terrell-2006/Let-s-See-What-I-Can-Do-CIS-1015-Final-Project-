# DNA Analyzer & Mutation Classifier

import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

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

# Mutation Classification 
def classify_mutation(seq1, seq2):
    mutations = []   
    
    for i in range(0, min(len(seq1), len(seq2)), 3):  
        codon1 = transcribe(seq1[i:i+3])   
        codon2 = transcribe(seq2[i:i+3])
        
        if codon1 != codon2 and len(codon1) == 3:   
            aa1 = codon_table.get[codon1, "?"]   
            aa2 = codon_table.get(codon2, "?")  
            
            if aa1 == aa2:   
                mtype = "Silent"   
            elif aa2 == "STOP":   
                mtype = "Nonsense"  
            else:
                mtype = "Missense"   
            
            mutations.append((i, codon1, codon2, aa1, aa2, mtype))  
    
    return mutations

# FASTA file reader 
def read_Fasta(filePath):
    sequence = ""   
    with open (filepath, "r") as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip()
    
    return sequence.upper()

# Main Analysis (Review Again) 
def analyze():
    seq1 = entry1.get().upper()
    
    if not validate_dna(seq1):
        output.delete("1.0", tk.END)
        output.insert(tk.END, "Invalid DNA sequence.\n")
        return
    
    rna = transcribe(seq1)
    protein = translate(rna)
    gc = gc_content(seq1)
    
    result = f"GC Content: {gc}%\nRNA: {rna}\nProtein: {protein}\n"
    
    seq2 = entry2.get().upper()
    if seq2:
        mutations = classify_mutations(seq1, seq2)
        result += "\nMutations:\n"
        for m in mutations:
            result += f"Pos {m[0]}: {m[1]}→{m[2]} | {m[3]}→{m[4]} ({m[5]})\n"
    
    output.delete("1.0", tk.END)
    output.insert(tk.END, result)


# Load File  
def load_file():  
    filepath = filedialog.askopenfilename()   
    seq = read_fasta(filepath)  
    entry1.delete(0, tk.END)
    entry1.insert(seq, 0)  


# Plot Button 
def plot():
    seq = entry1.get().upper   
    if validate_dna(seq):  
        plot_gc(seq)   

# GUI (Look Over For Bugs)
root = tk.Tk()
root.title("Genomic Analysis Tool")

tk.Label(root, text="DNA Sequence 1").pack()
entry1 = tk.Entry(root, width=50)
entry1.pack()

tk.Label(root, text="DNA Sequence 2 optional").pack()   
entry2 = tk.Entry(root, width=50)   
entry2.pack()   

tk.Button(root, text="Analyze", command=analyze).pack()  
tk.Button(root, text="Load FASTA File", command=load_file).pack()   
tk.Button(root, text="Plot GC Content", command=plot).pack()  

output = tk.Text(root, height=15, width=60)
output.pack()  

root.mainloop() 
