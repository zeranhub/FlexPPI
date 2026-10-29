# Academic Figure Skill Asset Confirmation (verified against existing native scripts)
# b/c -> manuscript Fig2 dual-metric color convention + LineTrend -> param inherit
# Colors encode metrics; line styles encode Full prediction / MD reference.
# Native matplotlib vectors; all residues retained, original validity masks applied.
from pathlib import Path
from hashlib import sha256
import json,shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
 'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'xtick.direction':'out','ytick.direction':'out',
 'xtick.major.width':.6,'ytick.major.width':.6,'legend.frameon':False,'axes.unicode_minus':False,
 'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':450,'savefig.bbox':None})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'revision_v107_soft_alpha_overlay'
INPUT=ROOT/'source_inputs'
RMSF_COLOR='#2166AC';HCHI_COLOR='#B87936'
RMSF_LIMIT=(0,2.4);HCHI_LIMIT=(0,4.5)
SASA_COLOR='#3A8D8C';SASA_LIMIT=(0,160)

def prepare():
    INPUT.mkdir(exist_ok=True)
    src_cam=OLD/'figure3/source_data/cam_all25.csv'
    ub_dir=ROOT.parent/'revision_v93_dual_case_entropy/ubiquitin_audit'
    original=[src_cam,ub_dir/'original_frozen_Full_A003R.npz',ub_dir/'original_frozen_labels_A003R.npz',OLD/'figure3/source_data/ubiquitin_72_position_profile_with_original_mask.csv']
    for p in original:shutil.copy2(p,INPUT/p.name)
    cam=pd.read_csv(INPUT/'cam_all25.csv',float_precision='round_trip').sort_values('CaMKI_native_position')
    frames=[]
    a=pd.DataFrame({'case':'CaMKI','sample_id':cam.sample_id,'pdb_id':'1MXE','residue_id':cam.CaMKI_native_position,
       'amino_acid':cam.amino_acid,'valid_rmsf':cam.valid_rmsf,'valid_hchi':cam.valid_chi_entropy,
       'MD_RMSF_A':cam.target_rmsf,'Full_RMSF_A':cam.pred_rmsf,'MD_Hchi_nats':cam.target_chi_entropy,'Full_Hchi_nats':cam.pred_chi_entropy})
    frames.append(a)
    with np.load(INPUT/'original_frozen_Full_A003R.npz',allow_pickle=False) as z,np.load(INPUT/'original_frozen_labels_A003R.npz',allow_pickle=False) as lab:
        assert lab['metric_names'].tolist()==['rmsf','chi_entropy','sasa']
        assert np.array_equal(z['target'],lab['values'],equal_nan=True) and np.array_equal(z['valid'],lab['mask'])
        ix=np.flatnonzero(z['chain']=='A');assert len(ix)==72
        a=pd.DataFrame({'case':'Ubiquitin','sample_id':'A003R','pdb_id':'2QHO','residue_id':z['author_residue_id'][ix].astype(int),
          'amino_acid':z['amino_acid'][ix],'valid_rmsf':z['valid'][ix,0],'valid_hchi':z['valid'][ix,1],
          'MD_RMSF_A':z['target'][ix,0].astype(float),'Full_RMSF_A':z['mean'][ix,0].astype(float),
          'MD_Hchi_nats':z['target'][ix,1].astype(float),'Full_Hchi_nats':z['mean'][ix,1].astype(float)})
        a=a.sort_values('residue_id').reset_index(drop=True)
    previous=pd.read_csv(INPUT/'ubiquitin_72_position_profile_with_original_mask.csv',float_precision='round_trip')
    assert np.array_equal(a.residue_id.to_numpy(),previous.residue_id.to_numpy())
    assert np.array_equal(a.valid_hchi.to_numpy(),previous.valid_hchi.to_numpy())
    for col in ['MD_Hchi_nats','Full_Hchi_nats']:
        mask=a.valid_hchi.to_numpy();assert np.array_equal(a.loc[mask,col].to_numpy(),previous.loc[mask,col].to_numpy())
    frames.append(a)
    result=pd.concat(frames,ignore_index=True)
    assert len(result)==97 and not result.duplicated(['case','residue_id']).any()
    result.to_csv(ROOT/'case_profiles_all97_with_masks.csv',index=False)
    stats={}
    for case in ['CaMKI','Ubiquitin']:
        z=result[result.case.eq(case)]
        stats[case]={'positions':len(z),'metrics':{}}
        for metric,unit,mask_col,lim in [('RMSF','A','valid_rmsf',RMSF_LIMIT),('Hchi','nats','valid_hchi',HCHI_LIMIT)]:
            q=z.loc[z[mask_col]];x=q['MD_'+metric+'_'+unit].to_numpy();y=q['Full_'+metric+'_'+unit].to_numpy()
            assert np.isfinite(x).all() and np.isfinite(y).all()
            assert min(x.min(),y.min())>=lim[0] and max(x.max(),y.max())<=lim[1]
            stats[case]['metrics'][metric]={'valid_positions':len(q),'masked_residues':z.loc[~z[mask_col],'residue_id'].astype(int).tolist(),
                'MAE':float(np.mean(abs(x-y))),'MD_range':[float(x.min()),float(x.max())],'Full_range':[float(y.min()),float(y.max())]}
    assert stats['CaMKI']['metrics']['Hchi']['valid_positions']==21
    assert stats['Ubiquitin']['metrics']['Hchi']['valid_positions']==66
    p=result[result.case.eq('Ubiquitin')&result.amino_acid.eq('P')]
    assert p.residue_id.tolist()==[19,37,38] and p.valid_hchi.all() and p.MD_Hchi_nats.eq(0).all()
    manifest={'sources':[{'original':str(p),'frozen':str(INPUT/p.name),'sha256':sha256(p.read_bytes()).hexdigest()} for p in original],
       'statistics':stats,'protocol':{'CaMKI':'existing family-heldout fold1, seed7, epoch40, corresponding single Full checkpoint',
        'Ubiquitin':'existing chemical_identity_corrected_v1 heldout fold0, seed7, epoch40, corresponding single Full checkpoint'},
       'new_model_or_MD_run':False,'color_semantics':'blue RMSF / ochre Hchi; solid Full / dashed MD',
       'shared_limits':{'RMSF_A':RMSF_LIMIT,'Hchi_nats':HCHI_LIMIT},'invalid_Hchi':'raw values retained; both MD and Full curves have NaN gaps under original mask'}
    (ROOT/'case_profile_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    return result,manifest

def load_data():
    if not (ROOT/'case_profiles_all97_with_masks.csv').exists():return prepare()
    return pd.read_csv(ROOT/'case_profiles_all97_with_masks.csv',float_precision='round_trip'),json.loads((ROOT/'case_profile_manifest.json').read_text(encoding='utf-8'))

def method_handles():
    return [Line2D([],[],color='#262626',lw=1.2,ls='-',label='Full'),Line2D([],[],color='#62676B',lw=.95,ls=(0,(3,2)),label='MD')]

def load_sasa():
    # Both cases use exact starting complexes and PyMOL dot_solvent=1,
    # probe1.4A/dot_density4; no Hchi mask is applied to geometric SASA.
    cam_source=ROOT/'sasa_audit'/'CaMKI_all25_PyMOL_sidechain_SASA.csv'
    ub_source=INPUT/'all119_input_environment.csv'
    INPUT.mkdir(exist_ok=True)
    assert ub_source.exists(), 'Frozen solvent-exposure input is missing'
    cam=pd.read_csv(cam_source,float_precision='round_trip')
    ub=pd.read_csv(INPUT/ub_source.name,float_precision='round_trip');ub=ub[ub.chain.eq('A')].sort_values('input_residue_id')
    a=pd.DataFrame({'case':'CaMKI','residue_id':cam.CaMKI_native_position.astype(int),
       'sidechain_SASA_A2':cam.sidechain_SASA_A2,'whole_residue_SASA_A2':cam.whole_residue_SASA_A2})
    b=pd.DataFrame({'case':'Ubiquitin','residue_id':ub.input_residue_id.astype(int),
       'sidechain_SASA_A2':ub.sidechain_SASA_A2,'whole_residue_SASA_A2':ub.whole_residue_SASA_A2})
    result=pd.concat([a,b],ignore_index=True)
    profiles,_=load_data()
    result=result.merge(profiles[['case','residue_id','amino_acid']],on=['case','residue_id'],validate='one_to_one')
    assert len(result)==97 and np.isfinite(result.sidechain_SASA_A2).all()
    assert result.sidechain_SASA_A2.between(*SASA_LIMIT).all()
    result.to_csv(ROOT/'case_SASA_all97_starting_complex.csv',index=False)
    return result

def draw_sasa(fig,canvas_mm,box,case):
    z=load_sasa();z=z[z.case.eq(case)].sort_values('residue_id')
    W,H=canvas_mm;left,top,w,h=box
    a=fig.add_axes([left/W,1-(top+h)/H,w/W,h/H])
    a.bar(z.residue_id,z.sidechain_SASA_A2,width=.78,color=SASA_COLOR,alpha=.60,edgecolor='none')
    a.set(xlim=(z.residue_id.min()-.5,z.residue_id.max()+.5),ylim=SASA_LIMIT,
          ylabel='SASA (Å²)',xlabel='CaMKI residue' if case=='CaMKI' else 'Ubiquitin residue')
    a.set_xticks([294,300,306,312,318] if case=='CaMKI' else [1,20,40,60,72]);a.set_yticks([0,150])
    a.tick_params(length=2.5,pad=2,width=.6);a.xaxis.labelpad=2.5;a.yaxis.labelpad=2.5
    a.yaxis.label.set_color(mpl.colors.to_rgba(SASA_COLOR,.60))
    return a,z

def draw_case(fig,canvas_mm,box,case,xlabel=True):
    frame,manifest=load_data();z=frame[frame.case.eq(case)].sort_values('residue_id')
    W,H=canvas_mm;left,top,w,h=box
    a=fig.add_axes([left/W,1-(top+h)/H,w/W,h/H]);b=a.twinx()
    # Both zero-based axes have fixed identical scales across the two cases.
    a.set_ylim(*RMSF_LIMIT);b.set_ylim(*HCHI_LIMIT)
    a.set_yticks([0,1,2]);b.set_yticks([0,2,4])
    a.set_ylabel('RMSF (Å)',color=RMSF_COLOR,labelpad=2.5)
    b.set_ylabel('Hχ (nats)',color=HCHI_COLOR,labelpad=2.5)
    for axis,col in [(a,RMSF_COLOR),(b,HCHI_COLOR)]:
        axis.tick_params(length=2.5,pad=2,width=.6);axis.tick_params(axis='y',colors=col)
        axis.spines['top'].set_visible(False)
    a.spines['left'].set_color(RMSF_COLOR);a.spines['right'].set_visible(False)
    b.spines['right'].set_visible(True);b.spines['right'].set_color(HCHI_COLOR);b.spines['left'].set_visible(False)
    rows=[]
    for metric,unit,mask,col,axis in [('RMSF','A','valid_rmsf',RMSF_COLOR,a),('Hchi','nats','valid_hchi',HCHI_COLOR,b)]:
        for method,lw,style,alpha in [('MD',.95,(0,(3,2)),.76),('Full',1.2,'-',1.)]:
            field=method+'_'+metric+'_'+unit;y=z[field].where(z[mask]).to_numpy();x=z.residue_id.to_numpy()
            marker=dict(marker='o',ms=1.8,mfc='white' if method=='MD' else col,mec=col,mew=.4 if method=='MD' else 0.) if metric=='Hchi' else {}
            # Markers preserve valid Hchi singletons flanked by undefined Ala/Gly gaps.
            axis.plot(x,y,color=col,lw=lw,ls=style,alpha=alpha,zorder=2 if method=='MD' else 3,**marker)
            assert np.isfinite(y).sum()==z[mask].sum()
            rows.extend({'case':case,'metric':metric,'unit':unit,'series':method,'residue_id':int(i),'value':float(v)} for i,v in zip(x,y) if np.isfinite(v))
    ticks=[294,300,306,312,318] if case=='CaMKI' else [1,20,40,60,72]
    a.set_xlim(z.residue_id.min()-.5,z.residue_id.max()+.5);a.set_xticks(ticks)
    if xlabel:a.set_xlabel('CaMKI residue' if case=='CaMKI' else 'Ubiquitin residue',labelpad=2.5)
    else:a.tick_params(labelbottom=False)
    return a,b,pd.DataFrame(rows)

def standalone():
    W,H=183.,105.
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
    fig.legend(handles=method_handles(),loc='upper center',bbox_to_anchor=(.5,.995),ncol=2,handlelength=2.,columnspacing=1.4,frameon=False)
    rows=[]
    for label,case,left,title in [('a','CaMKI',14,'CaM–CaMKI (1MXE)'),('b','Ubiquitin',109,'Ubiquitin–UBA (2QHO)')]:
        fig.text((left-11)/W,1-11/H,label,fontsize=11,fontweight='bold',va='top')
        fig.text((left+30)/W,1-11/H,title,fontsize=10,ha='center',va='top')
        a,b,p=draw_case(fig,(W,H),(left,23,60,48),case,xlabel=False);rows.append(p)
        draw_sasa(fig,(W,H),(left,76,60,13),case)
    values=pd.concat(rows,ignore_index=True);assert len(values)==368
    values.to_csv(ROOT/'standalone_all368_plotted_values.csv',index=False)
    fig.canvas.draw();r=fig.canvas.get_renderer();outside=[]
    for t in fig.findobj(mpl.text.Text):
        if t.get_visible() and t.get_text():
            bb=t.get_window_extent(r)
            if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:outside.append(t.get_text())
    assert not outside,outside
    for ext in ['png','pdf','svg']:fig.savefig(ROOT/f'case_RMSF_Hchi_dual_axis.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
    plt.close(fig)
    (ROOT/'case_profiles_QA.json').write_text(json.dumps({'positions':97,'plotted_values':368,'outside_text':outside,
       'same_RMSF_and_Hchi_limits_across_cases':True,'SASA_positions':97,'SASA_definition':'starting-complex sidechain PyMOL1.4A dot_density4',
       'same_SASA_limits_across_cases':True,'invalid_Hchi_gaps_retained':True,'legal_Pro_zeros_retained':True,'visual_review':'pending'},indent=2),encoding='utf-8')

if __name__=='__main__':
    prepare();standalone();print('Two cases: blue RMSF / ochre Hchi, solid Full / dashed MD; 368 original valid values.')
