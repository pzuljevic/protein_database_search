import string
from random import *

min_char = 50
max_char = 100
allchar = list("ACDEFGHIKLMNPQRSTVWY") 

def gen_protein_seq():
  return "".join(choice(allchar) for x in range(randint(min_char, max_char)))

# Returns similar protein sequence, where similar means
# there will be at most p% different charcters
def gen_similar_seq(seq, p):
  n = (int) (p * len(seq)) 
  similar_seq = list(seq)
  # Generate random positions to change
  pos = []
  for i in range (0, n): 
    similar_seq[randint(0, len(seq) - 1)] = allchar[randint(0, 19)]
  return ''.join(similar_seq)

K = 5 
seeds = [gen_protein_seq() for i in range(0, K)]
diff_pct = [0.01, 0.10, 0.25] 

for seed in seeds:
  l = str(len(seed))
  print ">C" + l
  print seed
  for p in diff_pct: 
    print ">C" + l + "_" + str(p)
    print gen_similar_seq(seed, p) 

