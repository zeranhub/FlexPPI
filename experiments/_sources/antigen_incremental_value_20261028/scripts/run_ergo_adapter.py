from pathlib import Path
import runpy
R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
adapter=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\ergo\run_ergo.py')
module=runpy.run_path(str(adapter),run_name='reused_verified_inference_adapter')
module['run'](R/'baseline_inputs.csv',R/'baselines/ergo')
