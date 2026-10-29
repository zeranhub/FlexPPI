from pathlib import Path
import json,shutil,numpy as np
out=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v86_dense_evidence\pymol_context')
p=out/'inputs/1A06_raw_official.pdb';lines=p.read_text().splitlines();protein=[s for s in lines if s.startswith('ATOM  ')];(out/'inputs/1A06_observed_native.pdb').write_text('\n'.join(protein)+'\nEND\n')
src=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v85_1MXE\pymol\inputs\1MXE_CaM_native4_147_CaMKI294_318.pdb');shutil.copy2(src,out/'inputs/1MXE_model_input_native.pdb')