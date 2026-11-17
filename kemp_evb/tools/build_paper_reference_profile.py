"""Plot saved OpenMM data against the reported Q6 barrier height (not a Q6 curve)."""
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]


def svg():
    source=(ROOT/'viewers/parameter_scan.html').read_text()
    data=json.loads(re.search(r'<script id="data" type="application/json">(.*?)</script>',source,re.S).group(1));row=data['baseline'];paper=data['published_reference']
    xs=[v-row['total_product_offset_kcal_mol'] for v in row['x']];ys=row['y'];mean=paper['barrier_kcal_mol'];sd=paper['replica_sd_kcal_mol']
    xmin,xmax=min(xs),max(xs);ymin=min(v for v in ys if v is not None)-3;ymax=max(max(v for v in ys if v is not None),mean+sd)+3
    x=lambda v:90+760*(v-xmin)/(xmax-xmin);y=lambda v:420-310*(v-ymin)/(ymax-ymin)
    s='<svg xmlns="http://www.w3.org/2000/svg" width="920" height="540" viewBox="0 0 920 540"><rect width="920" height="540" fill="white"/><g font-family="sans-serif" fill="#22313a"><text x="460" y="28" text-anchor="middle" font-size="19">Kemp: short OpenMM profile and published Q6 barrier</text>'
    s+='<path d="M95 50H125" stroke="#777" stroke-width="2" stroke-dasharray="5 5"/><text x="140" y="55" font-size="13">Grey: paper Des27.7-IN barrier 13.2 ± 0.6 kcal/mol (replica SD; not a full profile)</text><path d="M95 73H125" stroke="#1769aa" stroke-width="2"/><text x="140" y="78" font-size="13">Blue: our short OpenMM profile; conditional barrier 9.57 kcal/mol, unvalidated</text>'
    s+=f'<rect id="paper-barrier-sd" x="90" y="{y(mean+sd):.3f}" width="760" height="{y(mean-sd)-y(mean+sd):.3f}" fill="#888" opacity="0.16"/><path id="paper-barrier-line" data-value="13.2" d="M90 {y(mean):.3f}H850" stroke="#777" stroke-width="2" stroke-dasharray="5 5"/>'
    s+='<path d="M90 110V420H850" fill="none" stroke="#345"/>'
    yticks=sorted(set([0]+[round(ymin+(ymax-ymin)*i/5,1) for i in range(6)]))
    for tick in yticks:
        if tick!=0 and abs(y(tick)-y(0))<20:continue
        s+=f'<text x="80" y="{y(tick)+5:.3f}" text-anchor="end" font-size="13">{tick:g}</text>'
    xticks=[xmin+(xmax-xmin)*i/5 for i in range(6)]
    if xmin<=0<=xmax:xticks=[v for v in xticks if abs(x(v)-x(0))>45]+[0]
    for tick in sorted(xticks):s+=f'<text x="{x(tick):.3f}" y="445" text-anchor="middle" font-size="13">{tick:.1f}</text>'
    s+=f'<path d="M90 {y(0):.3f}H850" stroke="#cdd6dc" stroke-dasharray="3 4"/>'
    points=[]
    for gx,gy in zip(xs,ys):
        if gy is None:
            if points:s+='<polyline points="'+' '.join(points)+'" fill="none" stroke="#1769aa" stroke-width="2"/>';points=[]
            continue
        points.append(f'{x(gx):.3f},{y(gy):.3f}');s+=f'<circle cx="{x(gx):.3f}" cy="{y(gy):.3f}" r="3" fill="#1769aa"/>'
    if points:s+='<polyline points="'+' '.join(points)+'" fill="none" stroke="#1769aa" stroke-width="2"/>'
    s+='<text x="460" y="480" text-anchor="middle">Offset-adjusted EVB energy gap (kcal/mol)</text><text transform="translate(22 270) rotate(-90)" text-anchor="middle">Free energy relative to reactant (kcal/mol)</text><text x="460" y="510" text-anchor="middle" font-size="13">Paper: Listov et al., Nature2025, Supplementary Table3. Different Hamiltonians and sampling.</text></g></svg>'
    return s


def build():
    s=svg()
    (ROOT/'assets'/'source_profile.svg').write_text(s)
    p=ROOT/'viewers/source_profile.html';text=p.read_text();text=re.sub(r'<svg\b.*?</svg>',lambda m:s,text,flags=re.S)
    explanation='<p><strong>Grey is the published barrier height, not the paper’s profile curve.</strong> The shaded band is replica SD. Blue is our short, unvalidated OpenMM result. Short sampling and changed electrostatics/boundaries can both explain disagreement; longer sampling alone does not guarantee13.2.</p><p><a href="parameter_scan.html">Change H12 and alpha and watch the blue candidate</a> · <a href="https://doi.org/10.1038/s41586-025-09136-2" target="_blank" rel="noopener">Paper and supporting information</a></p>'
    text=re.sub(r'<div id="paper-reference-explanation">.*?</div>', '', text, flags=re.S)
    if explanation in text:text=text.replace(explanation,'')
    text=text.replace('</body>','<div id="paper-reference-explanation">'+explanation+'</div></body>');p.write_text(text)


if __name__=='__main__':build()
