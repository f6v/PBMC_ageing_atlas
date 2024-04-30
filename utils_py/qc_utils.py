import pandas as pd

def filter_counts(row):
  if row['total_counts'] > row['total_counts_min'] and \
    row['total_counts'] < row['total_counts_max'] and \
    row['pct_counts_mt'] < row['pct_counts_mt_max']:
    return 'pass'
  else:
    return 'fail'

def get_qc_df(adata, qc_cutoffs):
  qc_df = pd.merge(
    adata.obs,
    qc_cutoffs,
    left_on='sample_name',
    right_on='sample_name',
    how='left',
    sort=False
  )
  qc_df['barcode'] = adata.obs.index
  qc_df['qc'] = qc_df.apply(filter_counts, axis=1)
  
  print(qc_df['qc'].value_counts())
  
  return qc_df

def format_qc_df(qc_df):
	result =  qc_df.groupby(['sample_name', 'qc']).size().reset_index(name='count')
	result = result.pivot(index='sample_name', columns='qc', values='count').reset_index()
	result['total'] = result['fail'] + result['pass']
	result['failed_pct'] = 100 * (result['fail'] / result['total'])

	return result