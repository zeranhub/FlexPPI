from pathlib import Path
import pandas as pd
r=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\imrex')
p=pd.read_csv(r.parent/'binding_panel_inputs.csv')
x=pd.DataFrame({'case_id':p.case_id,'cdr3':p.CDR3_beta,'antigen.epitope':p.peptide})
x.to_csv(r/'input_imrex.csv',sep=';',index=False)
print(x.to_string(index=False))
