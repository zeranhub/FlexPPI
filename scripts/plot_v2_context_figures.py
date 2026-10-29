# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) grouped context bars -> BarAblation/param inherit (the source asset is an
# ablation template; its embedded data and labels are not reused).
# (b) paired context response -> LineTrend/param inherit.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","Liberation Sans"],"font.size":8,"axes.titlesize":8,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":8,"figure.titlesize":9,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.6,"xtick.direction":"out","ytick.direction":"out","xtick.major.width":0.6,"ytick.major.width":0.6,"legend.frameon":False})
# Academic Figure Skill Nature/Cell/Science Color Palette — COPY VERBATIM
CATEGORICAL=["#2166AC","#B2182B","#1B7837","#F1A340","#762A83","#666666"]
CATEGORICAL_EXTENDED=["#2166AC","#B2182B","#1B7837","#F1A340","#762A83","#666666","#4393C3","#D6604D","#5AAE61","#B35806","#9970AB","#999999"]
DIVERGING=["#2166AC","#F7F7F7","#B2182B"]
SEQUENTIAL=["#F7FBFF","#6BAED6","#08306B"]
ACCENT_RED="#B2182B";GREY="#999999";BLACK="#222222"
# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({"pdf.fonttype":42,"svg.fonttype":"none","savefig.bbox":"tight","savefig.dpi":300})
def save_cns_figure(fig,filename):
    fig.savefig(f"{filename}.pdf",bbox_inches="tight",dpi=300);fig.savefig(f"{filename}.png",bbox_inches="tight",dpi=300)

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];MM=1/25.4

def main():
    a=ROOT/'artifacts/evaluation_v2/itrace_context';out=ROOT/'output/flexppi_v2_figures';out.mkdir(exist_ok=True)
    d=pd.read_csv(a/'context_summary.csv'); d=d[(d.statistic=='spearman') & d.chain_role.isin(['peptide','tcr_beta','tcr_alpha','mhc','b2m'])]
    roles=['peptide','tcr_beta','tcr_alpha','mhc','b2m'];labels=['Peptide','TCR beta','TCR alpha','MHC','beta2m'];conditions=['full','partner_removed','isolated'];cl=['Full complex','Partner removed','Isolated chain'];colors=[CATEGORICAL[1],CATEGORICAL[0],GREY]
    fig,ax=plt.subplots(2,1,figsize=(183*MM,128*MM),sharex=True)
    for row,metric,title,ylabel in [(0,'rmsf','RMSF context sensitivity','Spearman ρ'),(1,'chi_entropy',r'$H_\chi$ context sensitivity','Spearman ρ')]:
        width=.24;x=np.arange(len(roles))
        for j,c in enumerate(conditions):
            vals=[];lo=[];hi=[]
            for role in roles:
                z=d[(d.chain_role==role)&(d.context==c)&(d.metric==metric)]
                vals.append(float(z['mean'].iloc[0]));lo.append(float(z['ci_low'].iloc[0]));hi.append(float(z['ci_high'].iloc[0]))
            ax[row].bar(x+(j-1)*width,vals,width,color=colors[j],label=cl[j],edgecolor='none')
            ax[row].errorbar(x+(j-1)*width,vals,yerr=[np.array(vals)-lo,np.array(hi)-vals],fmt='none',ecolor=BLACK,capsize=2,lw=.6)
        ax[row].set_ylabel(ylabel);ax[row].set_title(title,loc='left');ax[row].text(-.08,1.06,'ab'[row],transform=ax[row].transAxes,fontsize=9,fontweight='bold')
    ax[1].set_xticks(np.arange(len(roles)),labels);ax[0].legend(bbox_to_anchor=(1.01,1),loc='upper left');fig.text(.5,.02,'Bars: equal-weight complex mean · error bars: 95% complex bootstrap CI · fixed-coordinate input ablation, not apo/holo MD',ha='center',fontsize=6.5)
    fig.subplots_adjust(hspace=.38,bottom=.18,right=.78);save_cns_figure(fig,str(out/'fig_v2_itrace_context_ablation'));fig.savefig(out/'fig_v2_itrace_context_ablation.svg',dpi=300,bbox_inches='tight');plt.close(fig)

if __name__=='__main__':main()
