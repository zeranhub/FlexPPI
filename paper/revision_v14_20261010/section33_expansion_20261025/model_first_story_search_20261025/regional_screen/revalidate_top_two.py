"""Reuse the existing trajectory verifier without importing torch or writing elsewhere."""
import ast
import runpy
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent
CORE=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
SOURCE=OUT.parents[1]/'verify_case_trajectories.py'
NS={'__name__':'regional_trajectory_audit','__file__':str(SOURCE),'ROOT':OUT,'CORE':CORE,'np':np}
constants=runpy.run_path(str(CORE/'training_payload/src/flexppi/constants.py'))
NS.update({key:constants[key] for key in ['CHI_ATOMS','THREE_TO_ONE']})
label_tree=ast.parse((CORE/'training_payload/src/flexppi/data/labels.py').read_text(encoding='utf-8'))
needed=['_kabsch_fit','_chainwise_rmsf','_dihedral_degrees','_entropy']
exec(compile(ast.Module(body=[node for node in label_tree.body if isinstance(node,ast.FunctionDef) and node.name in needed],type_ignores=[]),'canonical_label_functions','exec'),NS)
tree=ast.parse(SOURCE.read_text(encoding='utf-8'))
nodes=[]
for node in tree.body:
    if isinstance(node,(ast.Import,ast.ImportFrom)):
        if isinstance(node,ast.ImportFrom) and node.module and node.module.startswith('flexppi'):
            continue
        nodes.append(node)
    elif isinstance(node,ast.FunctionDef):
        nodes.append(node)
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(SOURCE),'exec'),NS)
for sample in ['A004F','A007U']:
    NS['main'](sample,True)
