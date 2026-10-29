from pathlib import Path
import json,pandas as pd,zipfile
base=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010'); f=base/'section33_expansion_20261025/model_first_story_search_20261025/revision_v81_absolute_profiles/figure3'; p=f/'source_data'
a=pd.read_csv(p/'cam_md_state_contrast_all_common_residues.csv');c=a[a.metric=='hchi'];b=pd.read_csv(p/'cam_md_state_contrast_all_residues_five_blocks.csv');b=b[b.metric=='hchi'];x=c.contact_a>0;y=c.contact_b>0;values=[]
for mask in [x&~y,~x&y,x&y,~x&~y]:
    z=b[b.mature_cam_number.isin(c.loc[mask,'mature_cam_number'])].groupby('block')[['md_a','md_b']].mean();values.extend(z.to_numpy().ravel())
z=b[b.mature_cam_number.isin(c.loc[c.amino_acid=='M','mature_cam_number'])][['md_a','md_b']].to_numpy().ravel()
qa=json.loads((f/'figure_QA.json').read_text(encoding='utf-8'));qa['five_window_ranges_checked']={'contact_groups':[float(min(values)),float(max(values))],'all_nine_Met':[float(z.min()),float(z.max())],'all_intervals_within_axes':True};assert min(values)>1.3 and max(values)<2.2 and z.min()>1 and z.max()<3.4
(f/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
archive=base/'paper_v81_20261025_图3双状态绝对值源文件.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as out:
    for path in sorted(f.rglob('*')):
        if path.is_file():out.write(path,path.relative_to(f.parent))
print(json.dumps(qa['five_window_ranges_checked']))