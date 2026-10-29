"""Run the installed author's PaperBanana PolishAgent on the authentic v98 PNG.

The existing local credential stays in process memory and is never logged.
The initial suggestion is checked before allowing the image-generation step,
so a failed/zero-quota text request does not trigger a second API request.
"""
from pathlib import Path
import os,sys,asyncio,json,hashlib,base64,shutil,importlib.util
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parent
TOOLS=Path(r'C:\Users\DELL\.cache\flexppi-figure-tools')
REPO=TOOLS/'PaperBanana'
RUN=ROOT/'polish_v98_01'
INPUT=ROOT.parent/'revision_v98_structure_fill/figure3/fig3_dual_case_entropy.png'
KEY=TOOLS/'credentials/gemini-session.key'
sha=lambda b:hashlib.sha256(b).hexdigest()

def dump(name,value):
    (RUN/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=str),encoding='utf-8')

async def run():
    from utils import generation_utils
    from utils.config import ExpConfig
    from agents.polish_agent import PolishAgent,PLOT_POLISH_AGENT_SYSTEM_PROMPT
    import yaml
    from PIL import Image
    guard.install_gemini_failfast(generation_utils)
    config=yaml.safe_load((REPO/'configs/model_config.yaml').read_text(encoding='utf-8')) or {}
    defaults=config.get('defaults',{})
    main_model=defaults.get('main_model_name') or 'gemini-3.1-pro-preview'
    image_model=defaults.get('image_gen_model_name') or 'gemini-3.1-flash-image-preview'
    cfg=ExpConfig(dataset_name='PaperBananaBench',task_name='plot',split_name='custom_v98',exp_mode='polish',retrieval_setting='auto',main_model_name=main_model,image_gen_model_name=image_model,work_dir=RUN)
    agent=PolishAgent(exp_config=cfg,system_prompt=PLOT_POLISH_AGENT_SYSTEM_PROMPT)
    style=(RUN/'style_guides/neurips2025_plot_style_guide.md').read_text(encoding='utf-8')
    b64=base64.b64encode(INPUT.read_bytes()).decode('ascii')
    errors=[]
    original_router=generation_utils.call_model_with_retry_async
    async def observed_router(*args,**kwargs):
        try:return await original_router(*args,**kwargs)
        except Exception as error:
            errors.append(error);raise
    generation_utils.call_model_with_retry_async=observed_router
    provenance={'source':'https://github.com/dwzhu-pku/PaperBanana','source_commit':(REPO/'SOURCE_COMMIT.txt').read_text().strip(),'entrypoint':'agents.polish_agent.PolishAgent._generate_suggestions + PolishAgent.process','workflow':'author plot-image polish: VLM critique -> image refinement','main_model':main_model,'image_model':image_model,'input':str(INPUT),'input_sha256':sha(INPUT.read_bytes()),'aspect_ratio':'4:5','image_size':'2k','started_at':datetime.now(timezone.utc).isoformat(),'credential_values_not_logged':True,'existing_manuscript_and_v98_unchanged':True,'scope':'visual design trial; numerical fidelity must be reviewed before scientific use'}
    dump('provenance.json',provenance)
    print('Stage: PaperBanana VLM suggestions',flush=True)
    suggestions=await agent._generate_suggestions(b64,style)
    if not suggestions or not suggestions.strip():
        if errors:raise errors[-1]
        raise RuntimeError('PaperBanana returned no usable suggestions; image step was not requested.')
    (RUN/'suggestions.txt').write_text(suggestions,encoding='utf-8')
    dump('status.json',{'stage':'suggestions_saved','image_generated':False})
    async def cached_suggestions(*args,**kwargs):return suggestions
    agent._generate_suggestions=cached_suggestions
    print('Stage: PaperBanana original PolishAgent image refinement',flush=True)
    data={'path_to_gt_image':'v98.png','content':'Use the supplied real statistical figure; improve composition while preserving every scientific datum.','additional_info':{'rounded_ratio':'4:5','figure_size':'2K','image_size':'2k'}}
    result=await agent.process(data)
    image_data=result.get('polished_plot_base64_jpg')
    if not image_data:
        if errors:raise errors[-1]
        raise RuntimeError('PaperBanana produced no refined image.')
    raw=base64.b64decode(image_data)
    from io import BytesIO
    image=Image.open(BytesIO(raw));image.load()
    image.save(RUN/'paperbanana_polished_v98.png')
    result.pop('polished_plot_base64_jpg',None)
    dump('result_metadata.json',result)
    dump('status.json',{'stage':'complete','image_generated':True,'image_size':list(image.size),'output':'paperbanana_polished_v98.png','sha256':sha((RUN/'paperbanana_polished_v98.png').read_bytes())})
    print('Saved candidate 1: '+str(RUN/'paperbanana_polished_v98.png')+f' ({image.width}x{image.height})',flush=True)

if __name__=='__main__':
    RUN.mkdir(parents=True,exist_ok=True)
    if (RUN/'status.json').exists():raise SystemExit('Refusing to overwrite a recorded run. Use a new run directory.')
    (RUN/'data/PaperBananaBench/plot').mkdir(parents=True,exist_ok=True)
    (RUN/'style_guides').mkdir(exist_ok=True)
    shutil.copy2(INPUT,RUN/'data/PaperBananaBench/plot/v98.png')
    style='''User target: Nature Communications / Genome Biology scientific multi-panel case-study figure.
The six panels are a authentic protein ribbons together, b complete 25-residue CaMKI RMSF for MD/Full/Intra,
c complete 72-position ubiquitin profile with 66 valid Hchi MD/Full pairs, d all six Gln SASA/Hchi pairs,
e the complete chi1/chi2/chi3/chi4 2x2 grid, and f four full 83-record correlation groups.
Find worthwhile layout/spacing/alignment improvements from the supplied actual figure.
Keep the two structures large, together, and proportional; preserve purple target and gray partner.
Keep every scientific panel, data mark, line coordinate, numeric value, axis limit and unit unchanged.
Preserve n=66/52/23/10, rho=.50/.45/.52/.13 and the weak chi4 case.
Preserve black MD, blue Full, brown Intra, teal key residues. Black text, white background, readable labels.
Prefer concise typography, balanced real data areas, shared legends, aligned axes, and coherent whitespace.
No tables, gray microtext, decorative cards, borders, local molecular insets, mechanisms, invented data or structures.
Q41/Q62 is the environmental comparison; Q41/Q2 is the low/high joint-prediction comparison.
Full/Intra are two models on the same bound structure, never apo/bound experiments.
The output is a visual design trial; it will be independently checked against the original data.
'''
    (RUN/'style_guides/neurips2025_plot_style_guide.md').write_text(style,encoding='utf-8')
    secret=KEY.read_text(encoding='utf-8-sig').strip()
    if not secret:raise SystemExit('The configured credential file is empty; no request made.')
    os.environ['GOOGLE_API_KEY']=secret
    spec=importlib.util.spec_from_file_location('paperbanana_local_guard',TOOLS/'run_paperbanana.py')
    guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
    sys.stdout=guard.RedactedStream(sys.stdout,secret);sys.stderr=guard.RedactedStream(sys.stderr,secret)
    sys.path.insert(0,str(REPO))
    try:asyncio.run(run())
    except Exception as error:
        detail=guard.failure_summary(error)
        detail.update(image_generated=False,existing_figure_unchanged=True)
        dump('failure.json',detail);dump('status.json',{'stage':'failed','image_generated':False,'reason':detail.get('reason'),'status':detail.get('status')})
        print('Run failed: '+json.dumps(detail,ensure_ascii=True),flush=True)
        raise SystemExit(2)
