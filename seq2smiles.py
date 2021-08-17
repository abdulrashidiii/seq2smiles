#!/usr/bin/env python

Rs = {
    'A': 'N[C@@H](C)C(=O)',
    'R': 'N[C@@H](CCCN:C(:N):N)C(=O)',
    'N': 'N[C@@H](CC(=O)N)C(=O)',
    'D': 'N[C@@H](CC(=O)[O-])C(=O)',
    'C': 'N[C@@H](CS%s)C(=O)',
    'Q': 'N[C@@H](CCC(=O)N)C(=O)',
    'E': 'N[C@@H](CCC(=O)[O-])C(=O)',
    'G': 'NCC(=O)',
    'H': 'N[C@@H](Cc1nc[nH]c1)C(=O)',
    'I': 'N[C@@H](C(C)CC)C(=O)',
    'L': 'N[C@@H](CC(C)C)C(=O)',
    'K': 'N[C@@H](CCCC[NH3+])C(=O)',
    'M': 'N[C@@H](CCSC)C(=O)',
    'F': 'N[C@@H](Cc1ccccc1)C(=O)',
    'P': 'N1[C@@H](CCC1)C(=O)',
    'S': 'N[C@@H](CO)C(=O)',
    'T': 'N[C@@H](C(O)C)C(=O)',
    'W': 'N[C@@H](Cc1c2ccccc2[nH]c1)C(=O)',
    'Y': 'N[C@@H](Cc1ccc(O)cc1)C(=O)',
    'V': 'N[C@@H](C(C)C)C(=O)',
    'a': 'N[C@H](C)C(=O)',
    'r': 'N[C@H](CCCN:C(:N):N)C(=O)',
    'n': 'N[C@H](C(=O)N)C(=O)',
    'd': 'N[C@H](C(=O)[O-])C(=O)',
    'c': 'N[C@H](CS%s)C(=O)',
    'q': 'N[C@H](CC(=O)N)C(=O)',
    'e': 'N[C@H](CC(=O)[O-])C(=O)',
    'h': 'N[C@H](Cc1nc[nH]c1)C(=O)',
    'i': 'N[C@H](C(C)CC)C(=O)',
    'l': 'N[C@H](CC(C)C)C(=O)',
    'k': 'N[C@H](CCCC[NH3+])C(=O)',
    'm': 'N[C@H](CCSC)C(=O)',
    'f': 'N[C@H](Cc1ccccc1)C(=O)',
    'p': 'N1[C@H](CCC1)C(=O)',
    's': 'N[C@H](CO)C(=O)',
    't': 'N[C@H](C(O)C)C(=O)',
    'w': 'N[C@H](Cc1c2ccccc2[nH]c1)C(=O)',
    'y': 'N[C@H](Cc1ccc(O)cc1)C(=O)',
    'v': 'N[C@H](C(C)C)C(=O)',
    'NHE': 'N',
    'N_term': '[NH3+]',
    'C_term': '[O-]',
    'Ac': 'CC(=O)',
    'NME': 'NC',
    'A(S)PRO': 'N1[C@@H](C[C@H](N)C1)C(=O)',
    'A(R)PRO': 'N1[C@@H](C[C@@H](N)C1)C(=O)',
    'STY': 'N[C@@H](Cc1ccc(OS(=O)(=O)O)cc1)C(=O)',
    'O': 'N1[C@@H](C[C@@H](O)C1)C(=O)',
    'Z': 'N1[C@@H](CCC1(=O))C(=O)',
    'GLA': 'N[C@@H](C(C(=O)[O-])(C(=O)[O-]))C(=O)'
    }

def where(seq, n):
  idx = []
  for i in range(len(seq)):
    if seq[i] == n:
      idx.append(i)

  return idx

def read_connect(C_num, disul_connect):
  str_format = [0] * C_num
  conns = disul_connect.split(',')
  for conn in conns:
    idx = 9
    conn = conn.split('-')
    for i in conn:
      str_format[int(i)-1] = idx
      idx = idx - 1

  return tuple(str_format)

def seq2smiles(seq, disul_connect=None):
  try:
    C_idx = where(seq, 'C')
    C_num = len(C_idx)

    if disul_connect:
      str_format = read_connect(C_num, disul_connect)
    if disul_connect == None:
      str_format = tuple([''] * C_num)

    if len(seq) > 2:
  
      first = seq[0]
      last = seq[-1]
      seq = seq[1:-1]
    
      if first in ['Ac', 'Z', 'P', 'O']:
        smiles = Rs[first]
      else:
        smiles = Rs['N_term'] + Rs[first][1:]
    
      for aa in seq:
        smiles = smiles + Rs[aa]
    
      if last in ['NHE', 'NME']:
        smiles = smiles + Rs[last]
      else:
        smiles = smiles + Rs[last] + Rs['C_term']
  
    elif len(seq) == 2:
      first = seq[0]
      last = seq[1]

      if first in ['P', 'O']:
        smiles = Rs[first]
      else:
        smiles = Rs['N_term'] + Rs[first][1:]

      smiles = smiles + Rs[last] + Rs['C_term']
  
    elif len(seq) == 1:
      first = seq[0]

      if first in ['P', 'O']:
        smiles = Rs[first] + Rs['C_term']
      else:
        smiles = Rs['N_term'] + Rs[first][1:] + Rs['C_term']

    return smiles % str_format

  except KeyError:
    print('Error. Sequence must be valid.')
    exit(1)
