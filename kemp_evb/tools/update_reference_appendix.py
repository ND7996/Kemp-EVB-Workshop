"""Refresh the cited browser/PDF appendix and workshop integrity manifest."""
import hashlib
import html
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
SLIDES=[
 ('The EVB Hamiltonian',[
  'H = [[VR, H12], [H12, VP + alpha]]. The same coordinates are evaluated under reactant and product force-field descriptions.',
  'Alpha shifts the product diagonal; H12 is the off-diagonal coupling. Both have energy units, and neither is the activation barrier.',
  'Diagonalization gives E-minus/E-plus and mixing weights. At equal diagonals the splitting is 2|H12|.',
  'The website calculator demonstrates this energy calculation. Lambda is a sampling weight; a Morse alpha is a separate inverse-length parameter.'
 ],'Warshel and Weiss1980: 10.1021/ja00540a008'),
 ('Native Q6: preparation to profile',[
  'Start from the authors Amber-prepared coordinates, topology and FEP mapping; retain state-specific charges/bonded terms and Q boundary conditions.',
  'Qdyn performs restrained relaxation, heating and equilibration. The source protocol uses Berendsen coupling, leap-frog integration and SCAAS/LRF.',
  'Published sampling: 30ns equilibration, then51 mapping windows x200ps, with30 replicas per case.',
  'QFEP connects windows and reconstructs the EVB profile from native energy records. Original Des27 .en files are absent from the inspected deposit.'
 ],'Listov2025: 10.1038/s41586-025-09136-2; Gutierrez-Rus2025: 10.1021/jacs.4c09428'),
 ('EVBOpenMM: actual adaptation',[
  'Import mapped R/P systems for6330 particles using EVBOpenMM, converting units, LJ coefficients and Morse energy origins.',
  'The current driver relaxes lambda_R=1, heats50/150/300K, equilibrates, then traverses21 trial windows with window equilibration and production.',
  'Our direct-electrostatics droplet, wall, fixed outer protein and Langevin dynamics differ from the Q6 boundary and dynamics setup.',
  'EVBOpenMM analysis reconstructs the saved profile. Exact Q6 equivalence requires matched forces, restraints, sampling and analysis; more fitting alone cannot establish it.'
 ],'EVBOpenMM source: Martin-Floor/EVBOpenMM, revision dedf301eb9b0021a4be200327cea53f5b7f034c6'),
 ('What has actually been reproduced?',[
  'Reported Q6 Des27.7-IN:13.2 +/-0.6 kcal/mol. Recomputed short OpenMM estimate:9.572629 kcal/mol; difference=-3.627371.',
  'Two shortened native Q6/QFEP replicas: 9.82 and 8.10; mean=8.96, SD=1.216, SEM=0.860 kcal/mol (n=2).',
  'All65 shortened native stages completed normally; recovered/sibling templates are documented. Full-duration production and the30-replica result remain incomplete.',
  'Full-enzyme EVBOpenMM finite-energy/force checks passed. Saved-data repeatability does not establish convergence or parameter transfer.'
 ],'Evidence: data/reproduction_audit.json; full native instructions in EVB_MHETase/Kemp_EVB/docs'),
 ('Parameter history and calibration',[
  'Hong2018 SI Table10 reports H12=225, alpha=-324 kcal/mol; reactive parameter Tables10-21 and aqueous-reference target21.2 are available.',
  'Risso2020 explicitly reuses that older Q5.10 calibration. The2025 Q6 archive updates Morse/LJ terms and supplies100.4/-110.',
  'Do not combine older coupling/offsets with the updated force field. The solution-reference barrier21.2 is not the enzyme barrier13.2.',
  'All papers, supporting information, archives and review limitations are linked in the website parameter-history section.'
 ],'Hong2018: 10.1038/s41467-018-06305-y; Risso2020: 10.1039/D0SC01935F'),
 ('Grey is not the paper profile',[
  'Grey horizontal line/band: published Q6 barrier13.2 +/-0.6 (replica SD). It is not a full paper profile.',
  'Blue: selected computed profile, changing with H12/alpha. Optional amber: archived parameters on our OpenMM samples.',
  'The binned maximum is not a validated molecular transition state. Shifting it to x=0 changes the display origin only.',
  'Target matching is an illustration of parameter sensitivity. The native barrier has not been independently reproduced and cannot be claimed from a matched teaching curve.'
 ],'Q6 method: Bauer et al.2018, 10.1016/j.softx.2017.12.001; full bibliography on the website')]


def append_pdf(path):
    import pymupdf as fitz
    pdf=fitz.open(path)
    if len(pdf)>10:pdf.delete_pages(10,len(pdf)-1)
    for i,(title,bullets,citation) in enumerate(SLIDES,1):
        page=pdf.new_page(width=960,height=540)
        page.draw_rect(page.rect,color=None,fill=(.97,.975,.965))
        page.insert_text((55,40),f'OPTIONAL DETAIL {i}/{len(SLIDES)}',fontsize=11,color=(.09,.42,.45))
        n=page.insert_textbox(fitz.Rect(55,55,905,115),title,fontsize=28,fontname='hebo',color=(.12,.18,.21))
        if n<0:raise ValueError('Title overflow: '+title)
        y=125
        for bullet in bullets:
            n=page.insert_textbox(fitz.Rect(60,y,900,y+76),'- '+bullet,fontsize=17,fontname='helv',lineheight=1.25,color=(.12,.18,.21))
            if n<0:raise ValueError('Body overflow: '+title)
            y+=80
        n=page.insert_textbox(fitz.Rect(55,460,905,520),citation,fontsize=11,fontname='helv',color=(.32,.38,.42))
        if n<0:raise ValueError('Citation overflow: '+title)
    content=pdf.tobytes(garbage=4,deflate=True);pdf.close();Path(path).write_bytes(content)


def refresh():
    p=ROOT/'presentation.html';s=p.read_text()
    s=re.sub(r'\s*<!-- BEGIN REFERENCE APPENDIX -->.*?<!-- END REFERENCE APPENDIX -->','',s,flags=re.S)
    extra='\n<!-- BEGIN REFERENCE APPENDIX -->\n'
    for title,bullets,cite in SLIDES:
        extra+='<section class="slide"><div class="kicker">Optional technical detail</div><h2>'+html.escape(title)+'</h2><ul>'+''.join('<li>'+html.escape(b)+'</li>' for b in bullets)+'</ul><p class="source">'+html.escape(cite)+'</p></section>\n'
    extra+='<!-- END REFERENCE APPENDIX -->\n'
    s=s.replace('  </main>',extra+'  </main>',1);p.write_text(s)
    append_pdf(ROOT/'Kemp_EVB_workshop.pdf')
    update_manifest()


def update_manifest():
    p=ROOT/'workshop_manifest.json';d=json.loads(p.read_text());d['slide_count']=16;d['core_slide_count']=10;d['optional_detail_slides']=6
    d['talk_slides']=[1,10];d.pop('practical_slides',None)
    d['scope']='Standalone Kemp workshop with EVB Hamiltonian, cited parameter history and explicit native-Q6/OpenMM reconstruction limits.'
    d['files']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='workshop_manifest.json' and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.pyo')}
    (ROOT/'workshop_manifest.json').write_text(json.dumps(d,indent=2)+'\n')


if __name__=='__main__':refresh()
