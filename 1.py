import random
class MySeq:
    def __init__(self,seq,seq_type='DNA'):
        self.seq = seq 
        self.seq_type = seq_type 
    def set_seq_biotype(self,bt):
        biotype = bt.upper()
        if biotype == 'DNA' or biotype == 'RNA' or biotype == 'PROTEIN':
            self.seq_type = biotype
        else:
            print('Wrong biotype')    
    def print_sequence(self):
        print('Sequence: ',self.seq)
    
    def get_seq_biotype(self):
        return self.seq_type 
    def show_info_seq(self):
        print('Seqeuence:' + self.seq + 'Type: ' + self.seq_type)
    def count_occurences(self,seq_search):
        return self.seq.count(seq_search)

        
class MyNumSeq(MySeq):
    def __init__(self,num_seq,seq_type = 'numeric'):
        super().__init__(num_seq,seq_type)
    
    def set_seq_biotype(self,st):
        seq_type = st.upper()
        if seq_type == 'DNA' or seq_type == 'RNA' or seq_type == 'PROTEIN':
            self.seq_type = seq_type
        elif seq_type == 'NUMERIC' or seq_type == 'NUM':
            self.seq_type = 'NUMERIC'
        else:
            print('Wrong seq_type')

list_of_NumSeqs = []
for i in range(100):
    list_of_NumSeqs.append(MyNumSeq(str(int(random.random()*1000000))))
                           
for i in range(len(list_of_NumSeqs)):
    list_of_NumSeqs[i].print_sequence()