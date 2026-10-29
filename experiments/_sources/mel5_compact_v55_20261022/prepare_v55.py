from pathlib import Path
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'mel5_composite_v54_20261022'
FIG=ROOT/'figure4';FIG.mkdir(exist_ok=True)
s=(OLD/'plot_body.py').read_text(encoding='utf-8')
s=s.replace('W,H=183,214','W,H=183,172')
s=s.replace("iface['delta_Hchi']=iface.hchi_full_b-iface.hchi_full_a\nrank=iface.sort_values(['delta_Hchi','chain_a','position_a'],kind='stable').reset_index(drop=True)\nassert rank.iloc[0].chain_a=='A' and int(rank.iloc[0].position_a)==72\n",'')
start=s.index('# (d) Same ranked23')
end=s.index('# (e,f) Existing three contrasts',start)
s=s[:start]+s[end:]
s=s.replace('# (e,f) Existing three contrasts','# (d) Two stacked metrics in one joint-response panel. Existing three contrasts')
s=s.replace("heading('e','Recognition-loop responses to TCR and peptide changes',146)","heading('d','Joint backbone and side-chain responses',103)")
s=s.replace('axes(26,155,145,13)','axes(26,112,145,13)').replace('axes(174,155,2,13)','axes(174,112,2,13)')
s=s.replace("heading('f','pMHC side-chain responses to TCR and peptide changes',182)\n",'')
s=s.replace('axes(26,191,145,13)','axes(26,143,145,13)').replace('axes(174,191,2,13)','axes(174,143,2,13)')
s=s.replace("rank.to_csv(FIG/'source_interface_all23_ranked.csv',index=False,encoding='utf-8-sig')\n",'')
s=s.replace('panels=6,all_user_selected_panels_retained=True','panels=4,removed_panel="v54 d: ranked Hchi bars",merged_panels="v54 e and f -> v55 d upper and lower",all_requested_remaining_data_retained=True')
s=s.replace('ranked_Hchi=23,','ranked_Hchi=0,')
s=s.replace('fig4_mel5_composite','fig4_mel5_compact_joint')
(ROOT/'plot_body.py').write_text(s,encoding='utf-8')
prepare=(OLD/'prepare_figure.py').read_text(encoding='utf-8')
prepare=prepare.replace("# (d) Signed rank bars -> BarDistribution/plot_selfcorrection_math.py -> param inherit\n",'')
prepare=prepare.replace('# (e,f) Complete aligned delta maps','# (d upper/lower) Complete aligned delta maps')
prepare=prepare.replace("'BarDistribution/plot_selfcorrection_math.py',",'')
(ROOT/'prepare_figure.py').write_text(prepare,encoding='utf-8')

# Preserve the prior minimal OOXML builder and change only this revision's scope.
b=(OLD/'build_v54.py').read_text(encoding='utf-8')
head=b[:b.index('new={}\n')]
head=head.replace('six-panel composite','four-panel compact composite')
head=head.replace("SOURCE=BASE/'paper_v53_20261022_中文修订_MEL5残基接触与双动态联合图_三线表.docx'","SOURCE=BASE/'paper_v54_20261022_中文修订_MEL5六面板组合图4_三线表.docx'")
head=head.replace("OUT=BASE/'paper_v54_20261022_中文修订_MEL5六面板组合图4_三线表.docx'","OUT=BASE/'paper_v55_20261022_中文修订_MEL5精简联合图4_三线表.docx'")
head=head.replace('fig4_mel5_composite.png','fig4_mel5_compact_joint.png')
head=head.replace('d35736b5165c9e12525f414ffa2a2c966bfde17d01c7247310bff1834ebe9852','90b66a4029e22b0f67e52d3176839e547a8265690e5ff21d385ab7b5fecd1027')
new='''new={}
new[91]=old[91].replace('（图4d）','（图4d下部）')
new[93]=old[93].replace('（图4e）','（图4d上部）').replace('（图4f）','（图4d下部）')
new[95]='图4：MEL5体系中受体亲和成熟与肽替换的残基动态响应。（a）成熟复合物概览及Q72–β53/β54局部结构；虚线为Q72 OE1–β53 N和NE2–β54 O的原子间距。（b）四态全部10个肽位置的预测RMSF；P4为识别区域，P10为C端HLA锚定位点。（c）亲和成熟前后全部45个CDR位点的RMSF。（d）三种对照的联合动态响应：上部为CDR ΔRMSF，下部为pMHC界面ΔHχ；两项分别标示单位和色标。三行依次为4JFF−3HG1、4JFD−4JFF和4JFE−4JFF，蓝、红表示降低、升高，灰格表示Hχ无定义或氨基酸类型不同。界面为四态最近Cα距离≤8 Å的并集，共30位，三行Hχ有效数为23、23、22。TCR采用MEL5共同编号；各态均由同一模型对已有复合物结构预测，状态定义见方法。'
'''
tail=b[b.index('changes=[]\n'):]
tail=tail.replace('Figure 4 MEL5 six-panel structural and dynamic response composite','Figure 4 MEL5 compact four-panel joint response')
tail=tail.replace('Six panels: complex and Q72-beta53/54 atom-pair distances; four-state peptide RMSF; all45 CDR absolute RMSF; all23 valid interface Hchi changes; three-contrast complete CDR and pMHC interface heatmaps.','Four panels: complex and Q72-beta53/54 atom-pair distances; four-state peptide RMSF; all45 CDR absolute RMSF; unified upper CDR RMSF and lower pMHC Hchi heatmaps for three contrasts. No ranked bar panel.')
tail=tail.replace("assert all(s in case for s in ['图4b','图4c','图4d','图4e','图4f','G4A','I7A','5.85','10.58','2.94','3.12'])","assert all(s in case for s in ['图4b','图4c','图4d上部','图4d下部','G4A','I7A','5.85','10.58','2.94','3.12'])\nassert not any(s in alltext for s in ['图4e','图4f','（e、f）','p表示肽'])")
tail=tail.replace('prose_v54.json','prose_v55.json').replace('v54_content.txt','v55_content.txt')
(ROOT/'build_v55.py').write_text(head+new+tail,encoding='utf-8')
old_doc='paper_v54_20261022_中文修订_MEL5六面板组合图4_三线表.docx'
new_doc='paper_v55_20261022_中文修订_MEL5精简联合图4_三线表.docx'
for name in ['export_v54_word.ps1','render_v54.py']:
    t=(OLD/name).read_text(encoding='utf-8').replace(old_doc,new_doc).replace('v54','v55')
    (ROOT/name.replace('v54','v55')).write_text(t,encoding='utf-8')
print('Prepared compact4panel figure, targeted DOCX builder and render wrappers.')
