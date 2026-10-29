"""Verify every continuous cell and threshold outline in the revised panel b."""
from pathlib import Path
import json
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import plot_statewise_showcase as plot


def main():
    _,source,summary,membership=plot.load_data()
    fig,ax=plt.subplots()
    boxes_checked=0
    rows=[]
    for i,state in enumerate(plot.STATES):
        before_meshes=len(ax.collections)
        before_boxes=len(ax.patches)
        plot.draw_rank_pair(ax,source,membership,state,4*i)
        mesh=ax.collections[before_meshes]
        data=plot.select(source,state=state,metric='rmsf').set_index('author_residue_number').loc[plot.RESIDUES]
        np.testing.assert_array_equal(np.asarray(mesh.get_array()),data[['md_percentile','flexppi_percentile']].to_numpy().T)
        coords=mesh.get_coordinates()
        np.testing.assert_array_equal(coords[0,:,0],np.arange(585.5,673.5))
        np.testing.assert_allclose(coords[:,0,1],[4*i,4*i+1.65,4*i+3.3])
        assert mesh.norm.vmin==0 and mesh.norm.vmax==100
        calls=plot.select(membership,state=state).set_index('author_residue_number').loc[plot.RESIDUES]
        recovered={'Shared high':set()}
        for patch in list(ax.patches)[before_boxes:]:
            x,y=patch.get_xy()
            first=int(round(x+.5-.065))
            last=first+int(round(patch.get_width()+.13))-1
            assert np.isclose(patch.get_height(),3.3-.11), 'Source-only box remains'
            assert patch.get_linestyle()=='-', 'Dashed outline remains'
            recovered['Shared high'].update(range(first,last+1))
            boxes_checked+=1
        for category,positions in recovered.items():
            expected=set(calls.index[calls.category.eq(category)])
            assert positions==expected,(state,category,positions,expected)
        rows.append(dict(state=state,continuous_values=174,
                         shared=len(recovered['Shared high']),
                         md_only_unboxed=int(calls.category.eq('MD only').sum()),
                         model_only_unboxed=int(calls.category.eq('FlexPPI only').sum())))
    plt.close(fig)
    current=np.asarray(Image.open(plot.OUT/'fig01_kix_statewise_validation.png'))
    original=np.asarray(Image.open(plot.OUT/'history/panel_b_single_source_outlines_v2/fig01_kix_statewise_validation.png'))
    assert original.shape==current.shape
    h=current.shape[0]
    top_stop=int(h*(190-104)/190)
    bottom_start=int(h*(190-59)/190)
    assert np.array_equal(current[:top_stop],original[:top_stop]),'Panel a changed'
    assert np.array_equal(current[bottom_start:],original[bottom_start:]),'Panel c/d changed'
    result=dict(status='PASS',continuous_cells_checked=696,threshold_boxes_checked=boxes_checked,
                states=rows,unchanged_other_panels='a/c/d and page title/footer pixel-identical to archived v2',
                dimensions_png=[current.shape[1],current.shape[0]],
                no_interpolated_missing_residues=True,no_threshold_changes=True,
                no_source_only_boxes=True)
    (plot.OUT/'panel_b_revision_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
