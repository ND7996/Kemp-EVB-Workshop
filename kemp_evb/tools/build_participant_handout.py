from pathlib import Path
import pymupdf as fitz
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'PARTICIPANT_HANDOUT.pdf'
FIG=ROOT/'assets/kemp_states.png'
W,H=960,540
INK=(.12,.18,.21); MUTED=(.32,.39,.42); ACCENT=(.09,.42,.45); BLUE=(.13,.39,.68); GREEN=(.17,.48,.34); LINE=(.84,.88,.87); PALE=(.96,.975,.97); WHITE=(1,1,1)
doc=fitz.open()
def base(title,section):
 p=doc.new_page(width=W,height=H); p.draw_rect(p.rect,color=None,fill=PALE); p.draw_rect(fitz.Rect(30,25,W-30,H-25),color=LINE,fill=WHITE,width=1)
 p.insert_text((62,58),section.upper(),fontsize=10,fontname='hebo',color=ACCENT); p.insert_text((62,98),title,fontsize=29,fontname='hebo',color=INK)
 p.draw_line(fitz.Point(62,500),fitz.Point(W-62,500),color=LINE,width=.8)
 p.insert_text((62,519),'Kemp EVB workshop | Participant guide',fontsize=9,fontname='helv',color=MUTED)
 p.insert_text((W-92,519),f'{len(doc)}',fontsize=9,fontname='helv',color=MUTED)
 return p
def card(p,box,label,title,body,color=ACCENT):
 p.draw_rect(box,color=LINE,fill=WHITE,width=1); p.draw_rect(fitz.Rect(box.x0+9,box.y0+5,box.x1-9,box.y0+9),color=color,fill=color)
 p.insert_text((box.x0+17,box.y0+33),label.upper(),fontsize=10,fontname='hebo',color=color)
 p.insert_textbox(fitz.Rect(box.x0+17,box.y0+45,box.x1-17,box.y0+78),title,fontsize=16,fontname='hebo',color=INK)
 p.insert_textbox(fitz.Rect(box.x0+17,box.y0+88,box.x1-17,box.y1-13),body,fontsize=12,fontname='helv',color=INK,lineheight=1.25)

# Page 1: reaction
p=base('Kemp elimination: what changes?','Reaction at a glance')
p.insert_textbox(fitz.Rect(62,115,W-124,155),'Asp162 accepts a proton from 5-nitrobenzisoxazole as the substrate ring opens to form the cyanophenoxide product.',fontsize=14,fontname='helv',color=INK)
p.insert_image(fitz.Rect(90,155,W-90,425),filename=str(FIG),keep_proportion=True)
p.insert_textbox(fitz.Rect(62,432,924,488),'Atom labels: C3 is the substrate carbon that loses Ht; Ht is the transferred hydrogen; O1 and N2 are the ring oxygen and nitrogen. Asp162 is enzyme aspartate residue 162; OD1 is one of its side-chain oxygens and accepts Ht.',fontsize=11,fontname='helv',color=INK,lineheight=1.25)

# Page 2: the EVB idea and browser activities
p=base('Two bonding descriptions, one reaction','Understanding EVB')
boxes=[fitz.Rect(62,132,340,322),fitz.Rect(354,132,632,322),fitz.Rect(646,132,924,322)]
card(p,boxes[0],'Reactant state (R)','Before','The substrate still has its C-H bond and closed N-O-containing ring. Asp162 is a carboxylate ready to accept the proton.',BLUE)
card(p,boxes[1],'Product state (P)','After','Asp162 has accepted the proton. The N-O bond is broken, the nitrile has formed, and the ring-opened phenoxide product remains.',GREEN)
card(p,boxes[2],'EVB','How the states mix','EVB combines energy descriptions for R and P. The energy offset shifts one state relative to the other; coupling controls how strongly they mix.',ACCENT)
p.insert_text((62,365),'Explore the saved examples in this order',fontsize=18,fontname='hebo',color=INK)
activities=[('1','EVB fundamentals','Change coupling and offset one at a time.'),('2','Saved reaction profile','Find the reactant, high-energy region, and product.'),('3','Parameter comparison','Compare parameter pairs and their product energies.')]
for i,(n,t,b) in enumerate(activities):
 y=405+i*29; p.insert_text((68,y),n,fontsize=12,fontname='hebo',color=ACCENT); p.insert_text((94,y),t,fontsize=12,fontname='hebo',color=INK); p.insert_text((285,y),b,fontsize=11,fontname='helv',color=MUTED)
p.insert_textbox(fitz.Rect(62,478,924,496),'Paper Figure 3e shows barrier estimates; the OpenMM viewer shows a free-energy profile. Different boundaries, electrostatics, Morse references, and sampling mean these are not a direct comparison.',fontsize=9,fontname='helv',color=MUTED)

# Page 3: real calculation stages
p=base('From a structure to a free-energy profile','Full calculation workflow')
boxes=[fitz.Rect(62,132,340,358),fitz.Rect(354,132,632,358),fitz.Rect(646,132,924,358)]
card(p,boxes[0],'1 · Prepare','Build matched R and P','Use the same atoms in the same order. Assign state-specific bonds, charges, force-field terms, and EVB settings. Check the atom map, total charge, proton position, solvent, and ions before running.',BLUE)
card(p,boxes[1],'2 · Relax','Stabilize both endpoints','Minimize and equilibrate R and P separately. Inspect geometry, contacts, and energy stability. Relaxation prepares each endpoint; it does not calculate the reaction barrier.',ACCENT)
card(p,boxes[2],'3 · FEP/profile','Sample the lambda path','Sample a set of lambda windows between the state descriptions. Check neighboring-window overlap and convergence, compare independent replicas, then analyze the barrier and reaction free energy.',GREEN)
p.draw_rect(fitz.Rect(62,382,924,470),color=LINE,fill=PALE,width=1)
p.insert_text((80,409),'What this workshop contains',fontsize=15,fontname='hebo',color=INK)
p.insert_textbox(fitz.Rect(80,418,905,462),'This repository includes saved Kemp data and browser exercises for learning the method. It does not include the complete live R/P builder, endpoint-relaxation inputs, or lambda-window run scripts, so the three simulation stages above are explained but not launched here.',fontsize=12,fontname='helv',color=INK,lineheight=1.25)

doc.set_metadata({'title':'Kemp EVB participant guide','author':'Kemp EVB workshop'})
doc.save(OUT,garbage=4,deflate=True)
print(OUT)
