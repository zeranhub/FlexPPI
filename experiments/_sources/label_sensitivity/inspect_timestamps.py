import json
from pathlib import Path
from collections import Counter
root=Path(__file__).resolve().parent
rows=json.loads((root/'all_selected_inputs_before_window_definition.json').read_text(encoding='utf-8'))['records']
items=[item for row in rows for item in row['trajectories']]
print('records',len(rows),'trajectories',len(items))
print('steps',Counter((row['dataset'],item['median_step_ps']) for row in rows for item in row['trajectories']))
print('shortest_ps',min(item['last_ps'] for item in items))
for row in rows:
    for item in row['trajectories']:
        if item['duplicate_time_indices'] or item['backwards_time_indices']:
            print(row['sample_id'],item['replica_index'],'duplicates',item['duplicate_time_indices'],'backwards',item['backwards_time_indices'])
