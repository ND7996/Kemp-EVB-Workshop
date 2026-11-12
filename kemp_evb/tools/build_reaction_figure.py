from pathlib import Path
import math
import pymupdf as fitz

ROOT = Path(__file__).resolve().parents[1]
PNG = ROOT / 'assets' / 'kemp_states.png'
PDF = ROOT / 'assets' / 'kemp_reaction_scheme.pdf'
W, H = 1800, 760
INK = (0.12, 0.18, 0.21)
MUTED = (0.32, 0.39, 0.42)
BLUE = (0.13, 0.39, 0.68)
GREEN = (0.17, 0.48, 0.34)
TEAL = (0.09, 0.42, 0.45)
GOLD = (0.84, 0.48, 0.08)
RED = (0.72, 0.24, 0.20)
PALE = (0.96, 0.975, 0.97)
LINE = (0.84, 0.88, 0.87)
WHITE = (1, 1, 1)

doc = fitz.open()
p = doc.new_page(width=W, height=H)
p.draw_rect(p.rect, color=None, fill=PALE)

def text(x, y, value, size=24, color=INK, bold=False):
    p.insert_text((x, y), value, fontsize=size, fontname='hebo' if bold else 'helv', color=color)

def centered(x, y, value, size=22, color=INK, bold=False):
    font = fitz.Font('hebo' if bold else 'helv')
    tw = font.text_length(value, fontsize=size)
    text(x-tw/2, y, value, size, color, bold)

def pill(x, y, w, h, label, fill, fg=WHITE, size=19):
    p.draw_rect(fitz.Rect(x,y,x+w,y+h), color=None, fill=fill)
    centered(x+w/2, y+h/2+size*.34, label, size, fg, True)

def atom(x, y, label, fill, radius=28, size=22):
    p.draw_circle(fitz.Point(x,y), radius, color=None, fill=fill)
    centered(x,y+size*.34,label,size,WHITE,True)

def line(a,b,color=INK,width=5,dashes=None):
    p.draw_line(fitz.Point(*a),fitz.Point(*b),color=color,width=width,dashes=dashes)

def atom_bond(a,b,ra,rb,color=INK,width=5):
    dx,dy=b[0]-a[0],b[1]-a[1]
    length=math.hypot(dx,dy)
    ux,uy=dx/length,dy/length
    line((a[0]+ux*ra,a[1]+uy*ra),(b[0]-ux*rb,b[1]-uy*rb),color,width)

def point_to_atom(point,center,radius,color=INK,width=5):
    dx,dy=center[0]-point[0],center[1]-point[1]
    length=math.hypot(dx,dy)
    ux,uy=dx/length,dy/length
    line(point,(center[0]-ux*radius,center[1]-uy*radius),color,width)

def hexagon(cx,cy,r, nitro_side='left'):
    pts=[(cx,cy-r),(cx+r*.866,cy-r*.5),(cx+r*.866,cy+r*.5),(cx,cy+r),(cx-r*.866,cy+r*.5),(cx-r*.866,cy-r*.5)]
    for i in range(6): line(pts[i],pts[(i+1)%6],INK,5)
    # three inner lines suggest aromatic bonding without crowding atom labels
    for a,b in [(pts[0],pts[1]),(pts[2],pts[3]),(pts[4],pts[5])]:
        mx=(a[0]+b[0])/2; my=(a[1]+b[1])/2
        # offset inward toward ring centre
        dx=cx-mx; dy=cy-my; L=(dx*dx+dy*dy)**.5; off=(dx/L*13,dy/L*13)
        line((a[0]+off[0],a[1]+off[1]),(b[0]+off[0],b[1]+off[1]),MUTED,2)
    # nitro substituent on a benzene vertex; position is schematic
    v=pts[4] if nitro_side=='left' else pts[1]
    endpoint=(v[0]-48 if nitro_side=='left' else v[0]+48, v[1]+18)
    line(v,endpoint,INK,4)
    text(endpoint[0]-56,endpoint[1]+8,'NO2',20,BLUE,True)
    return pts

# Title and explanatory subtitle
text(64,72,'Kemp elimination: proton transfer opens the ring',38,INK,True)
text(66,112,'Asp162 accepts a substrate proton as the N-O bond breaks and the product forms.',23,MUTED)

# State cards
left=fitz.Rect(48,145,690,635); right=fitz.Rect(1110,145,1752,635)
p.draw_rect(left,color=LINE,fill=WHITE,width=2); p.draw_rect(right,color=LINE,fill=WHITE,width=2)
p.draw_rect(fitz.Rect(48,145,690,202),color=None,fill=(0.91,0.95,0.98))
p.draw_rect(fitz.Rect(1110,145,1752,202),color=None,fill=(0.91,0.96,0.93))
text(76,183,'REACTANT STATE',22,BLUE,True); text(76,226,'5-nitrobenzisoxazole',22,INK,True)
text(1138,183,'PRODUCT STATE',22,GREEN,True); text(1138,226,'2-cyano-4-nitrophenoxide',22,INK,True)

# Reactant schematic: fused aromatic ring and intact isoxazole ring
cx,cy,r=300,385,104
v=hexagon(cx,cy,r,'left')
# fused side v1-v2 and the three outer atoms O1-N2-C3
O=(445,279); N=(535,360); C3=(490,469)
point_to_atom(v[1],O,25,INK,5); atom_bond(O,N,25,25,RED,7); atom_bond(N,C3,25,25,INK,5); point_to_atom(v[2],C3,25,INK,5)
atom(O[0],O[1],'O1',RED,25,20); atom(N[0],N[1],'N2',BLUE,25,20); atom(C3[0],C3[1],'C3',INK,25,20)
# Acidic C-H bond, highlighted in gold
Ht=(585,502); atom_bond(C3,Ht,25,24,GOLD,6); atom(Ht[0],Ht[1],'Ht',GOLD,24,19)
text(564,540,'proton',16,GOLD,True)
# Asp catalytic base
p.draw_rect(fitz.Rect(84,552,620,610),color=LINE,fill=(0.97,0.98,0.99),width=1.5)
text(101,578,'Asp162',19,INK,True); text(101,601,'carboxylate',16,MUTED)
atom(370,581,'OD1',TEAL,23,16); text(400,588,'accepts Ht',16,TEAL,True)

# Middle reaction arrow and one-step label
center_x=900
pill(center_x-104,250,208,42,'ONE STEP',TEAL,size=19)
line((720,378),(1080,378),TEAL,8)
# arrow head
line((1080,378),(1054,362),TEAL,8); line((1080,378),(1054,394),TEAL,8)
centered(center_x,428,'concerted',19,TEAL,True)
centered(center_x,458,'bond changes',17,MUTED)
# Show concurrent changes underneath the arrow
centered(center_x,520,'C3-Ht  ->  OD1-Ht',18,GOLD,True)
centered(center_x,550,'N2-O1  breaks',18,RED,True)
centered(center_x,580,'C3-N2  becomes nitrile',18,GREEN,True)

# Product schematic: opened ring, phenoxide and nitrile
cx,cy,r=1330,385,104
v=hexagon(cx,cy,r,'left')
# O1 remains attached to aromatic fusion; N2 now forms a nitrile at C3
Oprod=(1460,260); point_to_atom(v[1],Oprod,25,INK,5); atom(Oprod[0],Oprod[1],'O1',RED,25,20)
text(1490,267,'-',22,RED,True)
Cprod=(1465,474); point_to_atom(v[2],Cprod,25,INK,5); atom(Cprod[0],Cprod[1],'C3',INK,25,20)
Nprod=(1560,474)
atom_bond(Cprod,Nprod,25,25,GREEN,8)
atom_bond((Cprod[0],Cprod[1]-9),(Nprod[0],Nprod[1]-9),20,20,GREEN,3)
atom_bond((Cprod[0],Cprod[1]+9),(Nprod[0],Nprod[1]+9),20,20,GREEN,3)
atom(Nprod[0],Nprod[1],'N2',BLUE,25,20)
# Asp protonated product state
p.draw_rect(fitz.Rect(1144,552,1680,610),color=LINE,fill=(0.97,0.98,0.99),width=1.5)
text(1161,578,'Asp162',19,INK,True); text(1161,601,'now protonated',16,MUTED)
atom(1470,581,'OD1',TEAL,23,16); text(1500,588,'- Ht',17,GOLD,True)

# Footer clarification
text(64,672,'Same atoms; different bonds and proton placement.',19,INK,True)
text(64,704,'Schematic only, not to scale. The paper describes a concerted reaction path; no separate stable intermediate is shown.',17,MUTED)

doc.set_metadata({'title':'Kemp elimination reaction scheme','author':'Kemp EVB workshop'})
doc.save(PDF,garbage=4,deflate=True)
pix=p.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False)
pix.save(PNG)
print(PNG)
print(PDF)
