#!/usr/bin/env python3
# rosalind.py

#%% import-libraries
from collections import Counter, deque
from functools import lru_cache

#%% define-functions
def dna(s):
    c = Counter(s)
    return "{} {} {} {}".format(c['A'], c['C'], c['G'], c['T'])

def rna(s):
    return s.replace('T', 'U')

def revc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

@lru_cache
def fib(n, k):
    if n < 3:
        return 1
    if k < 1:
        k = 1
    return k * fib(n - 2, k) + fib(n - 1, k)

def gc(s):
    d = deque(s.split('>'))
    d.popleft()
    kv = (x.split('\n', 1) for x in d)
    kv = ((k, ''.join(v.split())) for k, v in kv)
    kp = ((k, (v.count('C') + v.count('G')) / len(v) * 100) for k, v in kv)
    return "{}\n{:.6f}".format(*max(kp, key=lambda x: x[1]))

def hamm(s1, s2):
    return sum(x != y for x, y in zip(s1, s2))

def iprb(k, m, n):
    num = 4 * k * ((k - 1) + 2 * (m + n)) + m * (3 * (m - 1) + 4 * n)
    den = 4 * (k + m + n) * (k + m + n - 1)
    return "{:.5f}".format(num / den)

def prot(s):
    codon_table = {
        "UUU": 'F', "CUU": 'L', "AUU": 'I', "GUU": 'V',
        "UUC": 'F', "CUC": 'L', "AUC": 'I', "GUC": 'V',
        "UUA": 'L', "CUA": 'L', "AUA": 'I', "GUA": 'V',
        "UUG": 'L', "CUG": 'L', "AUG": 'M', "GUG": 'V',
        "UCU": 'S', "CCU": 'P', "ACU": 'T', "GCU": 'A',
        "UCC": 'S', "CCC": 'P', "ACC": 'T', "GCC": 'A',
        "UCA": 'S', "CCA": 'P', "ACA": 'T', "GCA": 'A',
        "UCG": 'S', "CCG": 'P', "ACG": 'T', "GCG": 'A',
        "UAU": 'Y', "CAU": 'H', "AAU": 'N', "GAU": 'D',
        "UAC": 'Y', "CAC": 'H', "AAC": 'N', "GAC": 'D',
        "UAA": '', "CAA": 'Q', "AAA": 'K', "GAA": 'E',
        "UAG": '', "CAG": 'Q', "AAG": 'K', "GAG": 'E',
        "UGU": 'C', "CGU": 'R', "AGU": 'S', "GGU": 'G',
        "UGC": 'C', "CGC": 'R', "AGC": 'S', "GGC": 'G',
        "UGA": '', "CGA": 'R', "AGA": 'R', "GGA": 'G',
        "UGG": 'W', "CGG": 'R', "AGG": 'R', "GGG": 'G'
    }
    codons = (s[i:i + 3] for i in range(0, len(s), 3))
    return ''.join(codon_table[codon] for codon in codons)

def subs(s1, s2):
    ii = range(len(s1) - len(s2))
    ns2 = len(s2)
    return ' '.join(str(i + 1) for i in ii if s1[i:i + ns2] == s2)

#%% tinkering...
def cons(s):
    d = deque(s.split('>'))
    d.popleft()
    kv = (x.split('\n', 1) for x in d)
    kv = ((k, ''.join(v.split())) for k, v in kv)
    return dict(kv)

#%% main
def main():
    # print(dna(""))
    # print(dna("A"))
    # print(dna("C"))
    # print(dna("G"))
    # print(dna("T"))
    # c = dna("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC")
    # print(c)
    
    # print(rna("GATGGAACTTGACTACGTAAATT"))
    
    # print(revc("AAAACCCGGT"))
    
    # print(fib(5, 3))
    
#     print(gc(""">Rosalind_6404
# CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
# TCCCACTAATAATTCTGAGG
# >Rosalind_5959
# CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
# ATATCCATTTGTCAGCAGACACGC
# >Rosalind_0808
# CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
# TGGGAACCTGCGGGCAGTAGGTGGAAT"""))
    
    # print(hamm("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"))
    
    # print(iprb(2, 2, 2))
    
    # print(prot("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
    
    # print(subs("GATATATGCATATACTT", "ATAT"))
    
    print(cons(""">Rosalind_1
ATCCAGCT
>Rosalind_2
GGGCAACT
>Rosalind_3
ATGGATCT
>Rosalind_4
AAGCAACC
>Rosalind_5
TTGGAACT
>Rosalind_6
ATGCCATT
>Rosalind_7
ATGGCACT"""))

if __name__ == "__main__":
    main()


