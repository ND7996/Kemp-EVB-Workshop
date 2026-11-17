from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Frame
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.colors import HexColor

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"PARTICIPANT_HANDOUT.pdf"
W,H=landscape((10*72,5.625*72))
ink=HexColor("#203039"); muted=HexColor("#53636b"); accent=HexColor("#176b73"); soft=HexColor("#eaf2f1")
styles=getSampleStyleSheet()
body=ParagraphStyle("Body",parent=styles["BodyText"],fontName="Helvetica",fontSize=13,leading=18,textColor=ink,spaceAfter=10)
title=ParagraphStyle("Title",parent=body,fontName="Helvetica-Bold",fontSize=27,leading=32,spaceAfter=18)
head=ParagraphStyle("Head",parent=body,fontName="Helvetica-Bold",fontSize=17,leading=21,textColor=accent,spaceAfter=8)
small=ParagraphStyle("Small",parent=body,fontSize=10,leading=14,textColor=muted)
def page(c,n,label,heading,items):
 c.setFillColor(HexColor("#f7f8f6"));c.rect(0,0,W,H,fill=1,stroke=0)
 c.setFillColor(accent);c.setFont("Helvetica-Bold",10);c.drawString(42,H-38,label.upper())
 p=Paragraph(heading,title);_,ph=p.wrap(W-84,H);p.drawOn(c,42,H-56-ph)
 frame=Frame(42,44,W-84,H-154,leftPadding=0,rightPadding=0,topPadding=6,bottomPadding=0,showBoundary=0)
 frame.addFromList(items,c)
 c.setStrokeColor(HexColor("#d8e1e1"));c.line(42,32,W-42,32)
 c.setFillColor(muted);c.setFont("Helvetica",9);c.drawString(42,19,"Kemp reaction and EVB hands-on");c.drawRightString(W-42,19,str(n))
 c.showPage()
c=canvas.Canvas(str(OUT),pagesize=(W,H))
c.setTitle("Kemp reaction and EVB participant guide")
page(c,1,"Start here","Three ways to model a reaction",[
 Paragraph("<b>Molecular dynamics (MD).</b> Moves atoms through time. In ordinary classical MD, bonds usually stay fixed.",body),
 Paragraph("<b>QM/MM MD.</b> Quantum mechanics describes the reacting atoms; a classical model describes the enzyme and solvent. It can represent bond changes, but costs more to run.",body),
 Paragraph("<b>EVB.</b> Couples classical force-field descriptions of reactant and product. With sampling, it can follow bond changes efficiently; parameters need suitable reference data.",body),
 Paragraph("The paper uses Q6 for EVB. This workshop uses OpenMM with EVBOpenMM for a short example. Both run EVB, with different setup and sampling.",small)
])
page(c,2,"Follow the chemistry","The Kemp reaction",[
 Paragraph("Asp162 in the enzyme accepts Ht from the substrate as the small ring opens.",body),
 Paragraph("<b>Atom labels:</b> C3 is the substrate carbon losing Ht; Ht is the transferred hydrogen; O1 and N2 are the ring atoms whose bond separates; Asp162 OD1 is the accepting oxygen.",body),
 Paragraph("Reactant and product are two bonding descriptions of the same atoms. In the 2 × 2 EVB table, the diagonal boxes hold the two state energies. α shifts the product energy; H₁₂ controls mixing in the other boxes. Neither is the reaction barrier.",body),
 Paragraph("A full calculation prepares matched states, relaxes them, samples FEP windows, then checks and analyzes the resulting profile.",body)
])
page(c,3,"Optional local run","Run the short example on Linux",[
 Paragraph("Install Git and Conda or Miniforge, then run from a terminal:",body),
 Paragraph("git clone https://github.com/ND7996/Kemp-EVB-Workshop.git<br/>cd Kemp-EVB-Workshop<br/>conda env create -f kemp_evb/calculation/environment.yml<br/>conda activate kemp-short<br/>python -m pip install --no-deps --no-build-isolation ./EVBOpenMM<br/>python kemp_evb/calculation/scripts/run_short_evb_demo.py --check<br/>python kemp_evb/calculation/scripts/run_short_evb_demo.py",small),
 Paragraph("The short run uses a CPU, five mapped atoms, and a prescribed geometry path. It writes an HTML plot and JSON summary. It is a software check, not molecular dynamics or a free-energy calculation.",body),
 Paragraph("Paper: Listov et al., Nature (2025), doi:10.1038/s41586-025-09136-2. Published simulation inputs: Zenodo 10.5281/zenodo.14563437.",small)
])
c.save()
print(OUT)
