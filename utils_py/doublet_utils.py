import vaeda
import pandas as pd
import numpy as np
from collections import Counter

def run_vaeda(adata, n_runs=10):
  all_vaeda_calls = pd.DataFrame()
  
  for i in range(0, n_runs):
    print('Run: %d' % i)
    adata = vaeda.vaeda(adata)
    all_vaeda_calls['vaeda_calls_%d' % i] = adata.obs['vaeda_calls']
    
    all_vaeda_calls['n_doublet_calls'] = all_vaeda_calls.apply(
      lambda row: Counter(row)['doublet'],
      axis=1
    )
    all_vaeda_calls['consensus_call'] = np.where(
      all_vaeda_calls['n_doublet_calls'] >= n_runs // 2,
      'doublet',
      'singlet'
    )
  
  return all_vaeda_calls