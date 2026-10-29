# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) True coordinate views -> render_biological_structures.py -> param inherit
# (b) All complete CDR/pMHC/peptide RMSF -> LineTrend visual logic -> param inherit
# (c) All same-AA pMHC interface Hchi -> PairedBoxScatter visual logic -> param inherit
# (d) All predefined contrast changes -> validated NPZ-derived CSV -> param inherit
# RULE: no native-run panels; retain every predefined biological contrast and position.
import os
os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from atlas_style import *
import pandas as pd
from matplotlib.lines import Line2D
ATLAS=ROOT.parent
DATA=ATLAS/'data_and_general_figures'
OUTROOT=ROOT/'biological_series'
OUTROOT.mkdir(exist_ok=True)
CAT=json.loads((ATLAS/'catalog/candidate_inventory.json').read_text(encoding='utf8'))
PALETTE=[FIRST,SECOND,'#2166AC','#1B7837','#B2182B','#666666']
MARKERS=['o','s','D','^','v','P']
CDRS=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
CONFIG={
 'C05':('PIK3CA H1047L recognition by two TCRs','interface',{'7l1d':'21LT2-2','7rrg':'0606T1-2'}),
 'C06':('MAGE-A10: c728 and its direct affinity-matured derivative c796','cdr',{'7pdw':'c728','7pbc':'c796'}),
 'C07':('AHIII12.2 recognition of three HLA single-site variants','cdr',{'1lp9':'WT','2j8u':'K66A','2jcc':'W167A','2uwe':'T163A'}),
 'C08':('1G4 recognition of NY-ESO-1 and its C9V anchor variant','peptide',{'2bnr':'C9','2bnq':'C9V'}),
 'C09':('N17.1.2 cross-recognition of NRAS Q61K and Q61R','cdr',{'8yiv':'Q61K','8yj2':'Q61R'}),
 'C10':('PRAME S2 recognition across HLA-A3 and HLA-A11','interface',{'8rym':'HLA-A3','8ryn':'HLA-A11'}),
 'C11':('KRAS G12V: 4TCR2 and the engineered MH derivative','cdr',{'8wte':'4TCR2','8wul':'MH'}),
 'C12':('PRAME S2 and the CDR1α-engineered S2-198 derivative','cdr',{'8rym':'S2','8ryo':'S2-198'}),
 'C13':('MEL5 recognition of MART-1 decamer and nonamer peptides','cdr',{'4qok':'Decamer','6eqa':'Nonamer'}),
 'C14':('DMF5 and DMF4 recognition of MART-1 peptide-length variants','parallel',{'3qdg':'DMF5 10mer','3qdj':'DMF5 9mer','3qdm':'DMF4 10mer','3qeq':'DMF4 9mer'}),
 'C15':('MEL5 recognition of Melan-A, IMP2 and BST2 antigens','cdr',{'4qok':'Melan-A','7q99':'IMP2','7q9a':'BST2'}),
 'C16':('1G4 and its CDR2-engineered c49c50 derivative','cdr',{'2f54':'WT 1G4','2f53':'c49c50'}),
 'C17':('TCR156 S32M: complete-complex readouts in a catch-bond study','whole',{'9nmu':'WT','9nmx':'S32M'}),
}
LABELS={'mhc':'MHC','peptide':'Peptide','tcr_alpha':'TCR α','tcr_beta':'TCR β','b2m':'β2M'}

def pair_frame(cid,a,b):
    path=DATA/cid/'data'/f'{a}_{b}_all_chain_alignment.csv'
    d=pd.read_csv(path,float_precision='round_trip',dtype={'id_a':str,'id_b':str})
    assert not d[['role','alignment_column']].duplicated().any()
    return d,path

def select_rmsf(d,scope):
    if scope=='cdr':
        m=d.region_a.isin(CDRS)|d.region_b_projected.isin(CDRS)
    elif scope in ['interface','parallel']:
        m=d.role.isin(['mhc','peptide'])&d.interface_union
    elif scope=='peptide':m=d.role.eq('peptide')
    else:m=d.paired&d.same_aa
    return d[m].copy()

def select_hchi(d,scope):
    if scope=='whole':m=d.paired&d.same_aa
    elif scope=='peptide':m=d.role.eq('peptide')&d.paired&d.same_aa
    else:m=d.role.isin(['mhc','peptide'])&d.interface_union&d.paired&d.same_aa
    return d[m&d.hchi_valid].copy()

def label_interface(d):
    return [('H' if r.role=='mhc' else 'P')+str(r.id_a) for r in d.itertuples()]

def finite_ticking(ax):
    for ticks,setter,lim in [(ax.get_xticks(),ax.set_xticks,ax.get_xlim()),(ax.get_yticks(),ax.set_yticks,ax.get_ylim())]:
        setter([t for t in ticks if lim[0]<=t<=lim[1]])

def plot_cdr(ax,frames,states,colors,markers):
    # All additional pairs share the same reference TCR in the multi-state cases.
    ref=frames[0][2]
    d=select_rmsf(ref,'cdr').copy()
    d['plot_region']=[a if a in CDRS else b for a,b in zip(d.region_a,d.region_b_projected)]
    keys=list(zip(d.role,d.id_a.fillna('gap'),d.alignment_column))
    arrays={}
    a0,b0,_=frames[0]
    arrays[a0]=d.rmsf_full_a.to_numpy(float)
    arrays[b0]=d.rmsf_full_b.to_numpy(float)
    for pa,pb,full in frames[1:]:
        if pa!=a0:continue
        t=select_rmsf(full,'cdr').copy()
        assert list(zip(t.role,t.id_a.fillna('gap'),t.alignment_column))==keys
        if pb not in arrays:arrays[pb]=t.rmsf_full_b.to_numpy(float)
        assert np.allclose(t.rmsf_full_a,d.rmsf_full_a,equal_nan=True)
    # Curves preserve all alignment rows. Missing or invalid values remain NaN.
    x=np.zeros(len(d));ticks=[];ticklabels=[];offset=0;records=[]
    for region in CDRS:
        ix=np.where(d.plot_region.to_numpy()==region)[0]
        assert len(ix),region
        xr=np.arange(len(ix))+offset;x[ix]=xr
        ax.text(float(xr.mean()),1.05,region,transform=ax.get_xaxis_transform(),ha='center',fontsize=9)
        ax.plot([xr[0]-.2,xr[-1]+.2],[1.01,1.01],transform=ax.get_xaxis_transform(),color='#B5B5B5',lw=.6,clip_on=False)
        for j in [ix[0],ix[-1]]:
            ticks.append(x[j]);ticklabels.append(d.iloc[j].id_a if not d.iloc[j].gap_a else '–/'+str(d.iloc[j].id_b))
        for state in arrays:
            y=arrays[state][ix]
            index=states.index(state)
            face='white' if index==0 else colors[state]
            ax.plot(xr,y,color=colors[state],marker=markers[state],ms=2.6,mfc=face,mec=colors[state],lw=.95)
        offset=xr[-1]+3
    for j,row in enumerate(d.itertuples()):
        for state,values in arrays.items():
            value=values[j]
            records.append({'region':row.plot_region,'plot_x':float(x[j]),'reference_id':row.id_a,
                            'state':state,'value_RMSF_A':None if not np.isfinite(value) else float(value)})
    values=np.concatenate([a[np.isfinite(a)] for a in arrays.values()])
    ylim=max(.4,np.ceil(values.max()*1.23*5)/5)
    ax.set(xlim=(-.6,float(x.max())+.6),ylim=(0,ylim),ylabel='RMSF (Å)',xticks=ticks)
    ax.set_xticklabels(ticklabels)
    ax.tick_params(length=2.5,pad=3)
    finite_ticking(ax)
    return d,records

def plot_pair_strip(ax,d,metric,colors,states,markers,rotate=60):
    a,b=states
    ac=metric+'_full_a';bc=metric+'_full_b'
    values=d[[ac,bc]].to_numpy(float)
    assert np.isfinite(values).all() and (values>=0).all()
    x=np.arange(len(d),dtype=float)
    for i,row in enumerate(d.itertuples()):
        va=float(getattr(row,ac));vb=float(getattr(row,bc))
        ax.plot([i-.14,i+.14],[va,vb],color='#AAAAAA',lw=.65,zorder=1)
        ax.plot(i-.14,va,marker=markers[a],mfc='white',mec=colors[a],ms=3.2,ls='none',zorder=3)
        ax.plot(i+.14,vb,marker=markers[b],mfc=colors[b],mec=colors[b],ms=3,ls='none',zorder=3)
    labels=label_interface(d)
    ymax=max(.4,np.ceil(values.max()*1.20*5)/5)
    ax.set(xlim=(-.65,len(d)-.35),ylim=(-max(.08,.06*ymax),ymax),xticks=x,
           ylabel='RMSF (Å)' if metric=='rmsf' else r'$H_\chi$ (nats)')
    ax.set_xticklabels(labels,rotation=rotate,ha='right' if rotate else 'center',rotation_mode='anchor')
    ax.tick_params(length=2,pad=2)
    finite_ticking(ax)
    assert ax.get_ylim()[1]>=values.max()
    return [{'plot_x':i,'role':r.role,'reference_id':r.id_a,'first':float(getattr(r,ac)),
             'second':float(getattr(r,bc)),'label':labels[i]} for i,r in enumerate(d.itertuples())]

def plot_peptide(ax,d,colors,states,markers):
    a,b=states
    x=np.arange(len(d))
    for state,col,val,valid,marker,face in [(a,colors[a],d.rmsf_full_a,d.rmsf_valid_a,markers[a],'white'),
        (b,colors[b],d.rmsf_full_b,d.rmsf_valid_b,markers[b],colors[b])]:
        y=val.to_numpy(float);y[~valid.to_numpy(bool)]=np.nan
        ax.plot(x,y,color=col,lw=1,marker=marker,ms=3.2,mfc=face,mec=col)
    labels=[]
    for r in d.itertuples():
        rid=r.id_a if not r.gap_a else r.id_b
        labels.append(f'{r.aa_a}{rid}' if r.same_aa else f'{r.aa_a if not r.gap_a else "–"}/{r.aa_b if not r.gap_b else "–"}{rid}')
    ymax=float(np.nanmax(d[['rmsf_full_a','rmsf_full_b']].to_numpy(float)))
    ax.set(xlim=(-.5,len(d)-.5),ylim=(0,max(.4,np.ceil(ymax*1.2*5)/5)),ylabel='RMSF (Å)',xticks=x)
    ax.set_xticklabels(labels)
    ax.tick_params(length=2.5,pad=3)
    finite_ticking(ax)
    return [{'plot_x':int(i),'label':labels[i],'first':None if r.gap_a else r.rmsf_full_a,
             'second':None if r.gap_b else r.rmsf_full_b} for i,r in enumerate(d.itertuples())]

def plot_whole(ax,full,metric,colors,states,markers):
    common=full.paired&full.same_aa&full[metric+'_valid']
    positions=np.arange(len(full),dtype=float)
    ticks=[];labels=[];records=[];offset=0
    ymax=float(full.loc[common,[metric+'_full_a',metric+'_full_b']].to_numpy(float).max())
    for role in ['mhc','peptide','tcr_alpha','tcr_beta','b2m']:
        d=full[full.role==role]
        if not len(d):continue
        ix=np.where(full.role.to_numpy()==role)[0]
        x=np.arange(len(ix),dtype=float)+offset
        rolemask=common.to_numpy()[ix]
        for suffix,state in zip(['a','b'],states):
            y=d[metric+'_full_'+suffix].to_numpy(float)
            y[~rolemask]=np.nan
            ax.plot(x,y,color=colors[state],lw=.65,marker=markers[state],ms=1.4,
                    mfc='white' if suffix=='a' else colors[state],mec=colors[state],mew=.35)
            for j,r in enumerate(d.itertuples()):
                records.append({'state':state,'role':role,'plot_x':float(x[j]),'reference_id':r.id_a,
                                'same_AA_valid':bool(rolemask[j]),'value':None if not rolemask[j] else float(y[j])})
        ax.text(float(x.mean()),1.05,LABELS[role],transform=ax.get_xaxis_transform(),ha='center',fontsize=9,color=GRAY)
        ax.plot([x[0],x[-1]],[1.01,1.01],transform=ax.get_xaxis_transform(),color='#B5B5B5',lw=.6,clip_on=False)
        role_order=[r for r in ['mhc','peptide','tcr_alpha','tcr_beta','b2m'] if (full.role==r).any()]
        choices=[] if role=='peptide' else [0,len(ix)//2]
        if role==role_order[-1]:choices.append(len(ix)-1)
        for j in sorted(set(choices)):
            ticks.append(x[j]);labels.append(str(d.iloc[j].id_a))
        offset=x[-1]+8
    ax.set(xlim=(-3,offset-5),ylim=(-max(.03,.05*ymax),max(.4,np.ceil(ymax*1.2*5)/5)),xticks=ticks,
           ylabel='RMSF (Å)' if metric=='rmsf' else r'$H_\chi$ (nats)')
    ax.set_xticklabels(labels)
    ax.tick_params(length=2.5,pad=3)
    finite_ticking(ax)
    return records,int(common.sum())

def diagnostic(ax,frames,scope,metric,colors,unit):
    allvalues=[];records=[];summaries=[]
    for contrast_index,(a,b,full) in enumerate(frames):
        d=select_rmsf(full,scope) if metric=='rmsf' else select_hchi(full,scope)
        if metric=='rmsf':d=d[d.paired&d.rmsf_valid]
        md=d[metric+'_md_delta'].to_numpy(float);pred=d[metric+'_full_delta'].to_numpy(float)
        assert len(md) and np.isfinite(md).all() and np.isfinite(pred).all()
        allvalues.extend(md.tolist()+pred.tolist())
        summaries.append({'a':a,'b':b,'n_positions':len(d),
                          'Spearman_rho':float(pd.Series(md).rank().corr(pd.Series(pred).rank())),
                          'MAE':float(np.mean(np.abs(md-pred))),
                          'zero_change_MAE':float(np.mean(np.abs(md))),
                          'direction_agreement':float(np.mean(np.sign(md)==np.sign(pred)))})
        d=d.copy();d['contrast']=a+'>'+b;records.append(d)
        ax.scatter(md,pred,s=12 if scope!='whole' else 8,c=colors[b],marker=MARKERS[contrast_index+1],
                   alpha=.75 if scope!='whole' else .45,edgecolors='none',zorder=3)
    extent=max(.05,max(abs(float(v)) for v in allvalues)*1.18)
    ax.plot([-extent,extent],[-extent,extent],lw=.6,color='#BBBBBB',ls='--',zorder=0)
    ax.axhline(0,color='#D0D0D0',lw=.5,zorder=0);ax.axvline(0,color='#D0D0D0',lw=.5,zorder=0)
    ax.set(xlim=(-extent,extent),ylim=(-extent,extent),xlabel=f'Existing MD change ({unit})',ylabel=f'Model change ({unit})')
    finite_ticking(ax)
    return pd.concat(records,ignore_index=True),summaries

def make_case(case):
    cid=case['case_id'];title,scope,names=CONFIG[cid]
    out=OUTROOT/cid;out.mkdir(exist_ok=True)
    frames=[];sources=[]
    for a,b in case['contrasts']:
        f,path=pair_frame(cid,a,b);frames.append((a,b,f));sources.append(path)
    a,b,first=frames[0]
    states=[]
    for pa,pb,_ in frames:
        for s in [pa,pb]:
            if s not in states:states.append(s)
    colors={s:PALETTE[i] for i,s in enumerate(states)}
    markers={s:MARKERS[i] for i,s in enumerate(states)}
    multi=len(frames)>1
    height=205 if multi or scope=='whole' else 200
    fig=make_canvas(height)
    heading(fig,'a',title,.975)
    assets=[]
    struct=DATA/cid/'structures'
    structure_view_path=struct/'molecular_view.json'
    structure_view=json.loads(structure_view_path.read_text(encoding='utf8'))
    for path,rect,lab,col in [
        (struct/f'{a}_overview.png',[.055,.806,.27,.143],f'{names.get(a,a.upper())} · {a.upper()}',colors[a]),
        (struct/f'{b}_overview.png',[.355,.806,.27,.143],f'{names.get(b,b.upper())} · {b.upper()}',colors[b]),
        (struct/'all_state_peptide_superposition.png',[.675,.824,.27,.10],'Aligned peptides',GRAY),
    ]:
        _,rec=add_molecule(fig,path,rect,lab,col);assets.append(rec)
    fig.text(.68,.802,f"All {len(case['pdbs'])} available structures",fontsize=9,color=GRAY)
    # All defined model states have a stable colour and marker across the figure.
    handles=[Line2D([],[],color=colors[s],marker=markers[s],mfc='white' if i==0 else colors[s],ms=3.6,lw=1,label=names.get(s,s.upper())) for i,s in enumerate(states)]
    fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.963,.786),ncol=min(len(states),4),
               columnspacing=1.2,handletextpad=.4,handlelength=1.3,borderaxespad=0)
    bdesc={'cdr':'backbone RMSF across all six complete CDRs','interface':'RMSF across the complete natural pMHC interface',
           'peptide':'RMSF across the complete aligned peptide','parallel':'complete pMHC RMSF in both receptor comparisons',
           'whole':'RMSF across all 407 unchanged complex positions'}[scope]
    heading(fig,'b','Model: '+bdesc,.737)
    brecords=[];hrecords=[];bframes=[];hframes=[]
    if scope=='parallel':
        for j,(pa,pb,f) in enumerate(frames):
            d=select_rmsf(f,scope);d=d[d.paired&d.rmsf_valid].copy()
            bottom=[.646,.507][j]
            ax=fig.add_axes([.105,bottom,.85,.045 if j==0 else .065])
            rowrecords=plot_pair_strip(ax,d,'rmsf',colors,[pa,pb],markers)
            for rec in rowrecords:rec.update({'a':pa,'b':pb})
            brecords.extend(rowrecords);bframes.append(d)
            label='DMF5' if pa=='3qdg' else 'DMF4'
            fig.text(.106,bottom+(.058 if j==0 else .077),f'{label}: {pa.upper()} → {pb.upper()} · n={len(d)}',fontsize=9,color=GRAY)
        heading(fig,'c','Model: all same-AA pMHC interface side chains in both TCRs',.439)
        for j,(pa,pb,f) in enumerate(frames):
            d=select_hchi(f,scope);bottom=[.336,.227][j]
            ax=fig.add_axes([.105,bottom,.85,.056])
            rr=plot_pair_strip(ax,d,'hchi',colors,[pa,pb],markers,rotate=0)
            for rec in rr:rec.update({'a':pa,'b':pb})
            hrecords.extend(rr);hframes.append(d)
            label='DMF5' if pa=='3qdg' else 'DMF4'
            fig.text(.106,bottom+.067,f'{label}: {pa.upper()} → {pb.upper()} · n={len(d)}',fontsize=9,color=GRAY)
        dhead=.183;dbottom=.064;dheight=.075
    elif scope=='whole':
        ax=fig.add_axes([.105,.530,.85,.157])
        brecords,nr=plot_whole(ax,first,'rmsf',colors,[a,b],markers)
        assert nr==407
        bframes=[first[first.paired&first.same_aa&first.rmsf_valid].copy()]
        fig.text(.53,.485,'Common unchanged positions; gaps retain undefined or substituted sites',fontsize=9,color=GRAY,ha='center')
        heading(fig,'c','Model: side-chain entropy at all 343 unchanged χ-defined sites',.453)
        ax=fig.add_axes([.105,.269,.85,.139])
        hrecords,nh=plot_whole(ax,first,'hchi',colors,[a,b],markers)
        assert nh==343
        hframes=[select_hchi(first,scope)]
        fig.text(.53,.226,'All complex roles retained · reference input residue numbering',fontsize=9,color=GRAY,ha='center')
        dhead=.195;dbottom=.064;dheight=.091
    else:
        ax=fig.add_axes([.105,.548 if not multi else .583,.85,.142 if not multi else .104])
        if scope=='cdr':
            bd,brecords=plot_cdr(ax,frames,states,colors,markers)
            bframes=[select_rmsf(f,scope) for _,_,f in frames]
            note='Reference TCR residue numbering; gaps separate complete CDRs'
        elif scope=='peptide':
            bd=select_rmsf(first,scope);brecords=plot_peptide(ax,bd,colors,[a,b],markers);bframes=[bd]
            note='Aligned peptide position; slash indicates altered amino-acid chemistry'
        else:
            bd=select_rmsf(first,scope);bd=bd[bd.paired&bd.rmsf_valid]
            brecords=plot_pair_strip(ax,bd,'rmsf',colors,[a,b],markers);bframes=[bd]
            note='Interface sites in sequence order; H: HLA/MHC, P: peptide'
        note_y=.479 if scope=='interface' else .495 if not multi else .546
        fig.text(.53,note_y,note,fontsize=9,color=GRAY,ha='center')
        if multi:
            heading(fig,'c','Model: all same-AA pMHC interface side chains by contrast',.514)
            bottoms=[.430,.328,.226]
            for j,(pa,pb,f) in enumerate(frames):
                d=select_hchi(f,scope);bottom=bottoms[j]
                ax=fig.add_axes([.105,bottom,.85,.048])
                rr=plot_pair_strip(ax,d,'hchi',colors,[pa,pb],markers,rotate=0 if len(d)<=20 else 45)
                for rec in rr:rec.update({'a':pa,'b':pb})
                hrecords.extend(rr);hframes.append(d)
                fig.text(.106,bottom+.059,f'{names.get(pa,pa.upper())} → {names.get(pb,pb.upper())} · n={len(d)}',fontsize=9,color=GRAY)
            dhead=.183;dbottom=.064;dheight=.075
        else:
            heading(fig,'c','Model: entropy at every comparable site in the defined region',.448)
            d=select_hchi(first,scope);hframes=[d]
            ax=fig.add_axes([.105,.308,.85,.105])
            hrecords=plot_pair_strip(ax,d,'hchi',colors,[a,b],markers)
            label='Complete same-AA peptide' if scope=='peptide' else 'Complete same-AA pMHC interface'
            fig.text(.53,.243,f'{label} · n={len(d)} · H: HLA/MHC, P: peptide',fontsize=9,color=GRAY,ha='center')
            dhead=.211;dbottom=.064;dheight=.107
    heading(fig,'d','Existing MD: all positions and all predefined contrasts',dhead)
    if multi:
        dhandles=[]
        for j,(pa,pb,_) in enumerate(frames):
            if scope=='parallel':lab='DMF5: 10mer→9mer' if pa=='3qdg' else 'DMF4: 10mer→9mer'
            else:lab=names.get(pa,pa.upper())+'→'+names.get(pb,pb.upper())
            dhandles.append(Line2D([],[],ls='none',marker=MARKERS[j+1],color=colors[pb],ms=3.6,label=lab))
        fig.legend(handles=dhandles,loc='upper right',bbox_to_anchor=(.963,dhead-.011),ncol=len(frames),
                   columnspacing=1.1,handletextpad=.35,handlelength=.8,borderaxespad=0)
    diagnostics={}
    for metric,rect,unit in [('rmsf',[.115,dbottom,.335,dheight],'Å'),('hchi',[.615,dbottom,.335,dheight],'nats')]:
        ax=fig.add_axes(rect)
        dd,ds=diagnostic(ax,frames,scope,metric,colors,unit)
        if multi:ax.set_ylabel('Model change\n('+unit+')',fontsize=9)
        ax.set_title(('Backbone' if metric=='rmsf' else 'Side chains')+f' · n={len(dd)}',fontsize=9,pad=1 if multi else 4)
        dd.to_csv(out/f'panel_d_{metric}_all_contrasts.csv',index=False)
        diagnostics[metric]=ds
    pd.DataFrame(brecords).to_csv(out/'panel_b_plotted_values.csv',index=False)
    pd.DataFrame(hrecords).to_csv(out/'panel_c_plotted_values.csv',index=False)
    for j,d in enumerate(bframes):d.to_csv(out/f'panel_b_scope_contrast{j+1}.csv',index=False)
    for j,d in enumerate(hframes):d.to_csv(out/f'panel_c_scope_contrast{j+1}.csv',index=False)
    nr=[int((d.paired&d.rmsf_valid).sum()) for d in bframes]
    nh=[int(len(d)) for d in hframes]
    nrep=10 if cid=='C17' else 3
    assert all(int(n)==nrep for _,_,f in frames for n in f.rmsf_nrep_a.dropna().unique())
    full_source_counts=[]
    for pa,pb,f in frames:
        br=select_rmsf(f,scope);hc=select_hchi(f,scope)
        full_source_counts.append({'a':pa,'b':pb,'backbone_scope_rows':len(br),
            'paired_valid_RMSF':int((br.paired&br.rmsf_valid).sum()),
            'Hchi_same_AA_valid':len(hc),'all_aligned_rows':len(f),
            'all_paired_rows':int(f.paired.sum()),'all_gap_rows':int((~f.paired).sum()),
            'all_altered_AA_rows':int((f.paired&~f.same_aa).sum())})
    audit=export_audited(fig,out/f'{cid}_case_figure',extra={
        'case_id':cid,'biological_contrasts':case['contrasts'],'all_structure_PDBs':case['pdbs'],
        'quantitative_states':states,'RMSF_scope':scope,'RMSF_valid_counts_by_contrast':nr,
        'Hchi_valid_counts_by_contrast':nh,'source_scope_counts':full_source_counts,
        'source_CSV_SHA256':{str(p):sha256(p) for p in sources},
        'source_audits':{str(p):sha256(p) for p in (DATA/cid/'data').glob('*source_audit.json')},
        'molecular_assets':assets,'diagnostics':diagnostics,'state_colors':colors,'state_markers':markers,
        'structure_view_manifest':{'path':str(structure_view_path),'SHA256':sha256(structure_view_path),
                                   'all_states':structure_view['all_states'],'colors':structure_view['semantic_colors'][:len(structure_view['all_states'])]},
        'all_predefined_contrasts_retained':True,'no_quality_based_selection':True,
        'no_smoothing_or_clipping':True,'primary_curves':'Absolute model predictions only; MD only in panel d',
        'gaps_and_undefined_chi':'Kept in complete CSV; NaN curve gaps or predefined valid-pair omission; never zero-filled',
        'n_MD_replicas_per_state':nrep,'new_training_inference_MD':False,
        'model_checkpoint_SHA256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',
    })
    color_names=['金橙','紫','蓝','绿','红','灰','浅蓝','棕','浅紫']
    overlay_mapping='、'.join(pdb.upper()+'（'+color_names[i]+'）' for i,pdb in enumerate(structure_view['all_states']))
    quantitative_names='、'.join(names.get(s,s.upper())+'（'+s.upper()+'）' for s in states)
    extras=[s for s in structure_view['all_states'] if s not in states]
    extras_text=('；'+ '、'.join(s.upper() for s in extras)+'仅出现在结构叠合中，不属于本图列明的定量对照') if extras else ''
    scope_label=('六个CDR' if scope=='cdr' else '共同未变全复合物' if scope=='whole' else '自然pMHC界面' if scope in ['interface','parallel'] else '肽段')
    counts_text='；'.join(names.get(pa,pa.upper())+'→'+names.get(pb,pb.upper())+f'：RMSF {nr[j]}位、Hχ {nh[j]}位' for j,(pa,pb,_) in enumerate(frames))
    role_text=('本例完整输入包括MHC、肽、TCRα和TCRβ，各角色按输入残基编号排列；343个侧链覆盖整个复合物，并非仅pMHC界面。' if scope=='whole' else 'H/P分别表示MHC/肽；CDR曲线间断处为不同环区段。' if scope=='cdr' else 'H/P分别表示MHC/肽。')
    caption=(f"{cid} {case['name']}。a，真实既有起始结构按MHC配准；左、中为{a.upper()}和{b.upper()}，右侧肽叠合的颜色依次对应{overlay_mapping}。定量状态为{quantitative_names}{extras_text}。b，完整预设{scope_label}的绝对模型RMSF。c，同氨基酸且χ角均有定义的全部{('肽段' if scope=='peptide' else '全复合物' if scope=='whole' else '自然pMHC界面')}侧链模型Hχ。各对照完整位数为{counts_text}。{role_text}d，在相同预设范围内比较模型变化与既有MD变化，统一为后态减前态；每个点为一个残基配对，虚线为数值一致。MD数值取各状态{nrep}次独立轨迹的等权均值；图中n指残基配对数，不能视作独立实验重复数。")
    (out/'caption_zh.txt').write_text(caption,encoding='utf8')
    stats={'case_id':cid,'status':'code_checks_passed_visual_review_pending','model':'Existing joint Full fold0; trained on ATLAS1866 + DynaRepo564; ITrace0',
           'model_n':1,'model_error_bars':None,'MD_independent_replicas_per_state':nrep,
           'MD_center':'Existing equal-weight replica mean','spread_displayed':None,
           'n_definition':'Aligned residue comparisons; not independent experimental samples or MD frames',
           'tests':None,'multiple_comparison_correction':None,'regression_fits':None,
           'descriptive_metrics_by_complete_contrast':diagnostics,
           'scope_counts':full_source_counts,'all_contrasts':case['contrasts'],
           'source_traceability':'Source CSVs and exact SHA256 in data_audit.json; panel_b/c plotted CSV and all panel_d observations exported',
           'limitations':case['model_cannot_support_or_limitations'],'PNG_SHA256':audit['PNG_SHA256']}
    (out/'statistics_and_reproducibility.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
    (out/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
    plt.close(fig)
    print(json.dumps({'case_id':cid,'RMSF_counts':nr,'Hchi_counts':nh,'output':str(out),'PNG_SHA256':audit['PNG_SHA256']},ensure_ascii=False),flush=True)
    return {'case_id':cid,'png':str(out/f'{cid}_case_figure.png'),'pdf':str(out/f'{cid}_case_figure.pdf'),
            'svg':str(out/f'{cid}_case_figure.svg'),'audit':str(out/'data_audit.json'),'caption':str(out/'caption_zh.txt')}

if __name__=='__main__':
    requested=set(sys.argv[1:])
    manifest=[]
    target=ROOT/'biological_series_manifest.json'
    previous=json.loads(target.read_text(encoding='utf8')) if target.exists() else {'figures':[]}
    for case in CAT['cases']:
        if case['case_id'] in CONFIG and (not requested or case['case_id'] in requested):
            manifest.append(make_case(case))
            keep=[x for x in previous['figures'] if x['case_id'] not in {x['case_id'] for x in manifest}]
            target.write_text(json.dumps({'status':'code_checks_passed_visual_review_pending','figures':sorted(keep+manifest,key=lambda x:x['case_id'])},ensure_ascii=False,indent=2),encoding='utf8')
