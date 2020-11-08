#!/usr/bin/env python3

#%%
from collections import Counter, deque
from functools import lru_cache

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
    pass

if __name__ == "__main__":
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
    
    print(prot("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
