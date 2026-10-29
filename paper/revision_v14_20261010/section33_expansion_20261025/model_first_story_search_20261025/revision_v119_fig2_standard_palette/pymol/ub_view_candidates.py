"""Re-render the exact input complex with different rigid camera views only."""
from pathlib import Path
import os,json,math,hashlib
import numpy as np
from PIL import Image,ImageDraw
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'view_candidates';OUT.mkdir(exist_ok=True)
FOV=20.
def project(xyz,view,w,h):
    v=np.asarray(view);q=(np.asarray(xyz)-v[12:15])@v[:9].reshape(3,3)+v[9:12]
    hh=-v[11]*math.tan(math.radians(FOV/2));hw=hh*w/h
    return np.column_stack([(0.5+q[:,0]/(2*hw))*w,(0.5-q[:,1]/(2*hh))*h])
def fit(c,w,h):
    v=list(c.get_view());R=np.array(v[:9]).reshape(3,3)
    xyz=c.get_coords('all');center=xyz.mean(0);q=(xyz-center)@R
    shift=(q.min(0)+q.max(0))/2;center=center+R@shift;q=(xyz-center)@R
    half=max(np.max(abs(q[:,1])),np.max(abs(q[:,0]))/(w/h))+2.2
    dist=half/math.tan(math.radians(FOV/2))
    v=list(R.ravel())+[0.,0.,-float(dist)]+list(center)+[1.,float(dist+200),-20.]
    c.set_view(v);return v
def render(c,name,w,h):
    v=fit(c,w,h)
    c.viewport(w,h);c.png(str(OUT/(name+'.png')),width=w,height=h,dpi=450,ray=1,quiet=1)
    xyz=c.get_coords('all');px=project(xyz,v,w,h)
    assert ((px[:,0]>=0)&(px[:,0]<w)&(px[:,1]>=0)&(px[:,1]<h)).all()
    lm={}
    for n in [41,62]:
        z=c.get_coords(f'chain A and resi {n} and name CA')
        lm['Q'+str(n)]=project(z,v,w,h)[0].tolist()
    return {'png':name+'.png','view':v,'landmarks':lm,'width':w,'height':h}
def main():
    configs=[('old',[]),('y45',[('y',45)]),('ym45',[('y',-45)]),('y90',[('y',90)]),('ym90',[('y',-90)]),
             ('x45',[('x',45)]),('xm45',[('x',-45)]),('yx',[('y',60),('x',30)]),
             ('orient',None),('orient_y45',None),('orient_ym45',None),('orient_x30',None)]
    all_entries={}
    with pymol2.PyMOL() as p:
        c=p.cmd;c.load(str(ROOT/'ub_overview.pse'))
        original=c.get_coords('all').copy();old=c.get_view();c.set('max_threads',4)
        c.set('field_of_view',FOV);c.set('orthoscopic',1)
        for name,turns in configs:
            c.set_view(old)
            if turns is None:
                c.orient('chain A')
                if name=='orient_y45':c.turn('y',45)
                if name=='orient_ym45':c.turn('y',-45)
                if name=='orient_x30':c.turn('x',30)
            else:
                for axis,deg in turns:c.turn(axis,deg)
            all_entries[name]=render(c,name,800,700)
            assert np.array_equal(original,c.get_coords('all'))
    (OUT/'candidates.json').write_text(json.dumps(all_entries,indent=2),encoding='utf-8')
    sheet=Image.new('RGB',(1600,3*400),'white');d=ImageDraw.Draw(sheet)
    for i,(name,z) in enumerate(all_entries.items()):
        im=Image.open(OUT/z['png']).convert('RGBA');im.thumbnail((390,350))
        x=(i%4)*400;y=(i//4)*400
        sheet.paste(im,(x+(400-im.width)//2,y+28),im)
        d.text((x+12,y+8),name,fill='black')
    sheet.save(OUT/'all_views.png')
    print('12 actual camera variants rendered; atomic coordinates unchanged.',flush=True)
if __name__=='__main__':main()
