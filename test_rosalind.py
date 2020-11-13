#!/usr/bin/env python3
# test_rosalind.py

#%% import-libraries
import rosalind
import unittest

#%% define-classes
class TestRosalind(unittest.TestCase):
    
    def test_dna(self):
        f = rosalind.dna
        self.assertEqual(f(''), "0 0 0 0")
        self.assertEqual(f('A'), "1 0 0 0")
        self.assertEqual(f('C'), "0 1 0 0")
        self.assertEqual(f('G'), "0 0 1 0")
        self.assertEqual(f('T'), "0 0 0 1")
        s = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
        self.assertEqual(f(s), "20 12 17 21")
    
    def test_rna(self):
        f = rosalind.rna
        self.assertEqual(f(''), '')
        self.assertEqual(f('A'), 'A')
        self.assertEqual(f('C'), 'C')
        self.assertEqual(f('G'), 'G')
        self.assertEqual(f('T'), 'U')
        self.assertEqual(f("GATGGAACTTGACTACGTAAATT"), "GAUGGAACUUGACUACGUAAAUU")
    
    def test_revc(self):
        f = rosalind.revc
        self.assertEqual(f(''), '')
        self.assertEqual(f('A'), 'T')
        self.assertEqual(f('C'), 'G')
        self.assertEqual(f('G'), 'C')
        self.assertEqual(f('T'), 'A')
        self.assertEqual(f("AAAACCCGGT"), "ACCGGGTTTT")
    
    def test_fib(self):
        f = rosalind.fib
        self.assertEqual(f(1, 1), 1)
        self.assertEqual(f(2, 1), 1)
        self.assertEqual(f(3, 1), 2)
        self.assertEqual(f(4, 1), 3)
        self.assertEqual(f(5, 1), 5)
        self.assertEqual(f(1, 99), 1)
        self.assertEqual(f(2, 99), 1)
        self.assertEqual(f(3, 2), 3)
        self.assertEqual(f(3, 3), 4)
        self.assertEqual(f(5, 3), 19)
    
    def test_gc(self):
        f = rosalind.gc
        s = """>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT"""
        r = """Rosalind_0808
60.919540"""
        self.assertEqual(f(s), r)
    
    def test_hamm(self):
        f = rosalind.hamm
        self.assertEqual(f('', ''), 0)
        self.assertEqual(f('A', 'A'), 0)
        self.assertEqual(f('A', 'C'), 1)
        self.assertEqual(f("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"), 7)
    
    def test_iprb(self):
        f = rosalind.iprb
        self.assertEqual(f(2, 0, 0), "1.00000")
        self.assertEqual(f(0, 2, 0), "0.75000")
        self.assertEqual(f(0, 0, 2), "0.00000")
        self.assertEqual(f(2, 2, 2), "0.78333")
    
    def test_prot(self):
        f = rosalind.prot
        s = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
        self.assertEqual(f(s), "MAMAPRTEINSTRING")
    
    def test_subs(self):
        f = rosalind.subs
        self.assertEqual(f("GATATATGCATATACTT", "ATAT"), "2 4 10")
    
    def test_cons(self):
        f = rosalind.cons
        s = """>Rosalind_1
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
ATGGCACT"""
        r = """ATGCAACT
A: 5 1 0 0 5 5 0 0
C: 0 0 1 4 2 0 6 1
G: 1 1 6 3 0 1 0 0
T: 1 5 0 0 0 1 1 6"""
        self.assertEqual(f(s), r)
    
    def test_fibd(self):
        f = rosalind.fibd
        self.assertEqual(f(1, 99), 1)
        self.assertEqual(f(99, 1), 0)
        self.assertEqual(f(99, 2), 1)
        self.assertEqual(f(2, 3), 1)
        self.assertEqual(f(3, 3), 2)
        self.assertEqual(f(4, 3), 2)
        self.assertEqual(f(5, 3), 3)
        self.assertEqual(f(6, 3), 4)
        self.assertEqual(f(6, 4), 6)
        self.assertEqual(f(6, 5), 7)
    
    def test_grph(self):
        f = rosalind.grph
        s = """>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG"""
        r = """Rosalind_0498 Rosalind_2391
Rosalind_0498 Rosalind_0442
Rosalind_2391 Rosalind_2323"""
        self.assertEqual(f(s), r)

if __name__ == "__main__":
    unittest.main()

