def validate_dna(dna_seq):
    seqm = dna_seq.upper()
    valid = seqm.count('A') + seqm.count('T') + seqm.count('G') + seqm.count('C')
    if valid == len(seqm):
        return True
    else:
        return False

def frequency(seq):
    dic = {}
    for s in seq.upper():
        if s in dic:
            dic[s] += 1
        else:
            dic[s] = 1
    return dic

def gc_content(seq):
    gc_count = 0 
    for s in seq:
        if s in 'GCgc':
            gc_count += 1
    return gc_count / len(seq)

def gc_content_subseq(dna_seq,k = 100):
    res = []
    for i in range(0,len(dna_seq)-k,k):
        subseq = dna_seq[i:i+k]
        gc = gc_content(subseq)
        res.append(gc)
    return res

def transcription(dna_seq):
    assert validate_dna(dna_seq), 'invalid dna sequence'
    return dna_seq.replace('T','U')

def reverse_complement(dna_seq):
    assert validate_dna(dna_seq), 'invalid dna sequence'
    comp = ""
    for c in dna_seq.upper():
        if c == 'A':
            comp += 'T'
        elif c == 'T':
            comp += 'A'
        elif c == 'G':
            comp += 'C'
        elif c == 'C':
            comp += 'G'
    return comp 

def translate_codon(cod):
    #make codon table
    codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}
    if cod in codontab:
        return codontab[cod]
    else:
        return None 

def translate_seq(dna_seq,ini_pos = 0):
    assert validate_dna(dna_seq), 'invalid dna sequence'
    seqm = dna_seq.upper()
    seq_aa = ''
    for pos in range(ini_pos,len(seqm)-2,3):
        cod = seqm[pos:pos+3]
        aa = translate_codon(cod)
        seq_aa += aa
    return seq_aa

def codon_usage(dna_seq,aa):
    assert validate_dna(dna_seq), 'invalid dna sequence'
    seqm = dna_seq.upper()
    dic = {}
    total = 0
    for i in range(0,len(seqm)-2,3):
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic:
                dic[cod] += 1
            else:
                dic[cod] = 1
            total += 1
    if total > 0:
        for k in dic:
            dic[k] /= total
    return dic

def reading_frames(dna_seq):
    assert validate_dna(dna_seq),'invalid dna sequence'
    res = []
    res.append(translate_seq(dna_seq,0))
    res.append(translate_seq(dna_seq,1))
    res.append(translate_seq(dna_seq,2))
    rc = reverse_complement(dna_seq)
    res.append(translate_seq(rc,0))
    res.append(translate_seq(rc,1))
    res.append(translate_seq(rc,2))
    return res

def all_proteins_rf(aa_seq):
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == '_':
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                    current_prot = []
        else:
            if aa == 'M':
                current_prot.append('')
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins 

def all_orfs(dna_seq):
    assert validate_dna(dna_seq), 'invalid dna sequence'
    rfs = reading_frames(dna_seq)
    res = []
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots:
            res.append(p)
    return res

def all_orfs_ord(dna_seq,minsize = 0):
    assert validate_dna(dna_seq), 'invalid dna sequence'
    rfs = reading_frames(dna_seq)
    res = []
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots:
            if len(p) > minsize:
                insert_prot_ord(p,res)
    return res

def insert_prot_ord(prot,list_prots):
    i = 0
    while i < len(list_prots) and len(prot) < len(list_prots[i]):
        i += 1
    list_prots.insert(i,prot)

class myseq:
    def __init__(self,seq,seq_type = 'DNA'):
        self.seq = seq.upper() 
        self.seq_type = seq_type 
    
    def __len__(self):
        return len(self.seq)

    def __getitem__(self,n):
        return self.seq[n]
    
    def __getslice__(self,i,j):
        return self.seq[i:j]
    
    def __str__(self):
        return self.seq
    
    def get_seq_biotype(self):
        return self.seq_type 
    
    def show_info_seq(self):
        print('Seqeuence:' + self.seq + 'Type: ' + self.seq_type)

    def alphabet(self):
        if (self.seq_type == 'DNA'):
            return 'ACGT'
        elif (self.seq_type == 'RNA'):
            return 'ACGU'
        elif (self.seq_type == 'PROTEIN'):
            return 'ACDEFGHIKLMNPQRSTVWY'
        else:
            return None
    
    def validate(self):
        alp = self.alphabet()
        res = True
        i = 0
        while i < len(self.seq) and res:
            if self.seq[i] not in alp:
                res = False
            else:
                i += 1
        return res

    def transcription(self):
        if (self.seq_type == 'DNA'):
            return self.seq.replace('T','U')
        else:
            return None
        
    def reverse_comp(self):
        if (self.seq_type != 'DNA'):
            return None 
        comp = ""
        for c in self.seq:
            if c == 'A':
                comp += 'T'
            elif c == 'T':
                comp += 'A'
            elif c == 'G':
                comp += 'C'
            elif c == 'C':
                comp += 'G'
        return myseq(comp,'DNA')
    
    def translate(self,iniPos = 0)
        if (self.seq_type != 'DNA'):
            return None
        seq_aa =''
        for pos in range(iniPos,len(self.seq)-2,3):
            cod = self.seq[pos:pos+3]
            aa = translate_codon(cod)
            seq_aa += translate_codon(cod)
        return myseq(seq_aa,'PROTEIN')
    

a = myseq('ATGCG',seq_type='RNA')
