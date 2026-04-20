# Validate DNA 
def validate_dna(seq):
    return all(base in "ATCG" for base in seq)

# Transcription 
def transcribe(seq):
    return seq.replace("T", "U")
