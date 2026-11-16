"""Build the short workshop slide PDF from the beginner hands-on outline."""
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph
from reportlab.lib.utils import ImageReader


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "Kemp_EVB_workshop.pdf"
PAGE_W, PAGE_H = 13.333 * inch, 7.5 * inch
INK = HexColor("#1f2e36")
MUTED = HexColor("#52626c")
ACCENT = HexColor("#176b73")
LINE = HexColor("#d7e0e3")
PAPER = HexColor("#f7f8f6")
SOFT = HexColor("#eaf2f1")

font_regular = Path(r"C:\Windows\Fonts\arial.ttf")
font_bold = Path(r"C:\Windows\Fonts\arialbd.ttf")
if font_regular.exists() and font_bold.exists():
    pdfmetrics.registerFont(TTFont("Workshop", str(font_regular)))
    pdfmetrics.registerFont(TTFont("Workshop-Bold", str(font_bold)))
    FONT, BOLD = "Workshop", "Workshop-Bold"
else:
    FONT, BOLD = "Helvetica", "Helvetica-Bold"

BODY = ParagraphStyle("Body", fontName=FONT, fontSize=16, leading=23, textColor=INK, spaceAfter=8)
LEAD = ParagraphStyle("Lead", parent=BODY, fontSize=22, leading=30)
SMALL = ParagraphStyle("Small", parent=BODY, fontSize=12, leading=17, textColor=MUTED)
CARD_BODY = ParagraphStyle("CardBody", parent=BODY, fontSize=14, leading=20)
CARD_TITLE = ParagraphStyle("CardTitle", parent=BODY, fontName=BOLD, fontSize=17, leading=22)
CENTER = ParagraphStyle("Center", parent=SMALL, alignment=TA_CENTER)


def para(c, text, x, top, width, style=BODY):
    p = Paragraph(text, style)
    _, height = p.wrap(width, PAGE_H)
    p.drawOn(c, x, top - height)
    return top - height


def begin(c, number, kicker, title):
    c.setFillColor(PAPER)
    c.rect(0, 0, PAGE_W, PAGE_H, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setStrokeColor(LINE)
    c.roundRect(34, 30, PAGE_W - 68, PAGE_H - 60, 8, fill=1, stroke=1)
    c.setFillColor(ACCENT)
    c.setFont(BOLD, 11)
    c.drawString(68, PAGE_H - 72, kicker.upper())
    title_style = ParagraphStyle("Title", fontName=BOLD, fontSize=34, leading=39, textColor=INK)
    y = para(c, title, 68, PAGE_H - 91, PAGE_W - 136, title_style)
    c.setStrokeColor(LINE)
    c.line(68, 51, PAGE_W - 68, 51)
    c.setFont(FONT, 10)
    c.setFillColor(MUTED)
    c.drawString(68, 35, "Kemp EVB workshop | Kemp reaction hands-on")
    c.drawRightString(PAGE_W - 68, 35, f"{number} / 11")
    return y - 18


def card(c, x, y, w, h, label, title, body, accent=None):
    c.setFillColor(colors.white)
    c.setStrokeColor(LINE)
    c.roundRect(x, y, w, h, 6, fill=1, stroke=1)
    if accent:
        c.setStrokeColor(accent)
        c.setLineWidth(4)
        c.line(x + 8, y + h - 3, x + w - 8, y + h - 3)
        c.setLineWidth(1)
    para(c, f'<font color="#176b73"><b>{label.upper()}</b></font>', x + 16, y + h - 17, w - 32, SMALL)
    para(c, f"<b>{title}</b>", x + 16, y + h - 52, w - 32, CARD_TITLE)
    para(c, body, x + 16, y + h - 84, w - 32, CARD_BODY)


def build():
    c = canvas.Canvas(str(OUTPUT), pagesize=(PAGE_W, PAGE_H), pageCompression=1)
    c.setTitle("Try EVB with the Kemp reaction")
    c.setAuthor("Kemp EVB workshop")

    # 1: Welcome
    y = begin(c, 1, "Browser-based hands-on", "Try EVB with the Kemp reaction")
    c.setFillColor(SOFT)
    c.roundRect(68, 100, 360, 50, 6, fill=1, stroke=0)
    c.setFillColor(ACCENT)
    c.setFont(BOLD, 14)
    c.drawString(86, 119, "Workshop website: nd7996.github.io/Kemp-EVB-Workshop")
    c.linkURL("https://nd7996.github.io/Kemp-EVB-Workshop/", (68, 100, 800, 150), relative=0)
    c.showPage()

    # 2: Reaction
    y = begin(c, 2, "1 · The reaction", "What changes?")
    y = para(c, "The substrate is the molecule the enzyme acts on. Asp162, part of the enzyme, accepts a proton (a hydrogen nucleus) from it. As the C–H bond breaks, the small ring opens and the product forms.", 68, y, 790, BODY)
    image = ImageReader(str(ROOT / "assets" / "kemp_states.png"))
    iw, ih = image.getSize()
    max_w, max_h = 810, 210
    scale = min(max_w / iw, max_h / ih)
    draw_w, draw_h = iw * scale, ih * scale
    c.drawImage(image, (PAGE_W - draw_w) / 2, 130, width=draw_w, height=draw_h, preserveAspectRatio=True, mask="auto")
    para(c, "Atom labels: C3 is the substrate carbon that loses Ht; Ht is the transferred hydrogen; O1 and N2 are the ring oxygen and nitrogen. Asp162 is enzyme aspartate residue 162; OD1 is one of its carboxylate oxygens and accepts Ht.", 68, 118, PAGE_W - 136, SMALL)
    c.showPage()

    # 3: Two states
    y = begin(c, 3, "2 · Before and after", "Two descriptions of the same atoms")
    box_y, box_h, gap = 175, 190, 22
    box_w = (PAGE_W - 136 - gap) / 2
    card(c, 68, box_y, box_w, box_h, "Reactant state", "Before", "The substrate holds the proton and its ring is closed. Asp162 is ready to accept it.", HexColor("#39769a"))
    card(c, 68 + box_w + gap, box_y, box_w, box_h, "Product state", "After", "Asp162 holds the proton. The substrate ring is open.", HexColor("#4c8768"))
    para(c, "EVB uses both descriptions to model the bond changes.", 68, 145, PAGE_W - 136, LEAD)
    c.showPage()

    # 4: EVB matrix in plain language
    y = begin(c, 4, "3 · EVB in one table", "How does the small matrix work?")
    para(c, "Think of the matrix as a 2 × 2 score table. It compares the reactant and product bond descriptions for the same atoms.", 68, y, PAGE_W - 136, LEAD)
    table_x, table_y, cell_w, cell_h = 92, 193, 196, 62
    c.setFillColor(MUTED); c.setFont(BOLD, 12)
    c.drawCentredString(table_x + cell_w / 2, table_y + cell_h * 2 + 15, "Reactant description")
    c.drawCentredString(table_x + cell_w * 1.5, table_y + cell_h * 2 + 15, "Product description")
    c.drawRightString(table_x - 10, table_y + cell_h * 1.5 - 4, "Reactant")
    c.drawRightString(table_x - 10, table_y + cell_h * .5 - 4, "Product")
    labels = [["Reactant energy", "H12: mixing"], ["H12: mixing", "Product energy + alpha"]]
    for row in range(2):
        for col in range(2):
            x = table_x + col * cell_w
            yy = table_y + (1 - row) * cell_h
            c.setFillColor(SOFT if row == col else colors.white)
            c.setStrokeColor(LINE); c.rect(x, yy, cell_w, cell_h, fill=1, stroke=1)
            c.setFillColor(INK); c.setFont(BOLD, 14)
            c.drawCentredString(x + cell_w / 2, yy + cell_h / 2 - 5, labels[row][col])
    para(c, "The diagonal boxes are the two state energies. Alpha shifts the product energy up or down. H12 in the other boxes controls how strongly the two descriptions mix.", 535, 354, 340, BODY)
    para(c, "At each geometry, the calculation combines the two energies and uses the lower mixed value. The table is a calculation tool, not two molecules. Alpha and H12 are not the reaction barrier.", 535, 244, 340, BODY)
    para(c, "Real values need suitable reference data and a validated model. The demonstration settings are not a fit.", 535, 165, 340, SMALL)
    c.showPage()

    # 5: Q6 and OpenMM
    y = begin(c, 5, "4 · Read the comparison", "These plots show different things")
    gap, box_y, box_h = 20, 185, 205
    box_w = (PAGE_W - 136 - gap) / 2
    card(c, 68, box_y, box_w, box_h, "Paper · Figure 3e", "Q6 activation barriers", "Compares experimental and Q6 barrier estimates for Des27 and Des27.7, including the in and out substrate conformations. It is not a free-energy profile trace.")
    card(c, 68 + box_w + gap, box_y, box_w, box_h, "This workshop · OpenMM", "Short-run free-energy profile", "The viewer plots saved data from 21 windows, 50 records per window, and 18 bins. Its 9.57 kcal/mol peak estimate is unvalidated; convergence has not been established.")
    c.linkURL("https://media.springernature.com/lw685/springer-static/image/art%3A10.1038%2Fs41586-025-09136-2/MediaObjects/41586_2025_9136_Fig3_HTML.png", (68, 189, 68 + box_w, 217), relative=0)
    c.linkURL("https://nd7996.github.io/Kemp-EVB-Workshop/viewers/source_profile.html", (68 + box_w + gap, 189, PAGE_W - 68, 217), relative=0)
    c.setFillColor(ACCENT)
    c.setFont(BOLD, 12)
    c.drawString(84, 198, "Open complete Figure 3")
    c.drawString(68 + box_w + gap + 16, 198, "Open the OpenMM profile viewer")
    para(c, "These plots report different quantities. Q6 used OPLS-AA, TIP3P, SCAAS solvent, and a local reaction field. Boundary treatment, electrostatics, Morse energy reference, and sampling differ in the OpenMM adaptation, so it is not a direct reproduction or engine-to-engine comparison.", 68, 155, PAGE_W - 136, BODY)
    para(c, "The Q6 paper reports activation-barrier estimates for its enzyme variants. A complete, converged OpenMM EVB calculation would produce a free-energy profile, activation free energy, and reaction free energy. The saved OpenMM peak is only an unvalidated teaching result.", 68, 100, PAGE_W - 136, SMALL)
    c.showPage()

    # 6: Full calculation workflow
    y = begin(c, 6, "5 · A full Kemp calculation", "Prepare, relax, then sample")
    gap, box_y, box_h = 14, 151, 240
    box_w = (PAGE_W - 136 - 2 * gap) / 3
    card(c, 68, box_y, box_w, box_h, "Prepare", "Build R and P", "Use the same atoms and atom order. Assign bonds, charges, force-field terms, and EVB settings. Check mapping, total charge, proton position, solvent, and ions.")
    card(c, 68 + box_w + gap, box_y, box_w, box_h, "Relax", "Settle both endpoints", "Minimize and equilibrate R and P separately. Check geometry and stability. This prepares the states; it does not calculate the barrier.")
    card(c, 68 + 2 * (box_w + gap), box_y, box_w, box_h, "FEP/profile", "Sample lambda windows", "Sample the path between R and P. Check window overlap and convergence, compare replicas, then analyze the profile, barrier, and reaction free energy.")
    para(c, "A reduced five-atom live demo is included with its mapped inputs and Slurm job. The separate full-enzyme backend contains the longer preparation-to-profile workflow.", 68, 125, PAGE_W - 136, SMALL)
    c.showPage()

    # 7: Core example
    y = begin(c, 7, "6 · Your turn", "Try changing the settings")
    para(c, "Open the small EVB example. Move one control at a time. Watch how the curve changes when you change the energy offset, then when you change the coupling.", 68, y, 780, LEAD)
    c.setFillColor(ACCENT)
    c.roundRect(68, 214, 310, 52, 6, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont(BOLD, 16)
    c.drawString(88, 233, "Open the EVB example")
    c.linkURL("https://nd7996.github.io/Kemp-EVB-Workshop/viewers/evb_core.html", (68, 214, 378, 266), relative=0)
    para(c, "The saved viewer needs no installation. A reduced five-atom live calculation is also included; see the cluster guide. It is not the full enzyme calculation.", 68, 184, 700, SMALL)
    c.showPage()

    # 8: The labeled profile
    y = begin(c, 8, "7 · The saved profile", "Kemp free-energy profile")
    image = ImageReader(str(ROOT / "assets" / "source_profile_slide.png"))
    c.drawImage(image, 68, 100, width=PAGE_W - 136, height=330, preserveAspectRatio=True, anchor="c", mask="auto")
    para(c, "This saved short-run profile is not validated.", 68, 90, PAGE_W - 136, SMALL)
    c.showPage()

    # 9: Parameter scan
    y = begin(c, 9, "8 · Compare settings", "One barrier does not choose one model")
    para(c, "Try the parameter comparison. Different offset/coupling pairs can give similar peak heights while giving different product-versus-reactant energies.", 68, y, 780, LEAD)
    c.setFillColor(ACCENT)
    c.roundRect(68, 214, 340, 52, 6, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont(BOLD, 16)
    c.drawString(88, 233, "Open the settings comparison")
    c.linkURL("https://nd7996.github.io/Kemp-EVB-Workshop/viewers/parameter_scan.html", (68, 214, 408, 266), relative=0)
    para(c, "These are teaching examples made from saved data. Matching a target is not a parameter fit or a validated prediction.", 68, 184, 740, SMALL)
    c.showPage()

    # 10: Take-home
    y = begin(c, 10, "9 · The point of the exercise", "What should you take away?")
    items = [
        "EVB connects reactant and product bonding descriptions.",
        "The offset and coupling shape the profile in different ways.",
        "Reference data and adequate sampling are needed before interpreting a profile physically.",
    ]
    for i, item in enumerate(items):
        yy = y - i * 47
        c.setFillColor(SOFT)
        c.circle(81, yy - 5, 7, fill=1, stroke=0)
        para(c, item, 102, yy + 8, 710, BODY)
    para(c, "Listov et al., Nature (2025). DOI: 10.1038/s41586-025-09136-2", 68, 140, 720, SMALL)
    c.linkURL("https://doi.org/10.1038/s41586-025-09136-2", (68, 48, 760, 70), relative=0)
    # 11: Cluster commands
    y = begin(c, 11, "10 · Run the live demo", "Submit a short calculation on Slurm")
    para(c, "Create the environment and verify the bundled inputs on the login node, from the workshop repository root:", 68, y, 800, BODY)
    command_style = ParagraphStyle("Command", parent=BODY, fontName="Courier", fontSize=13, leading=19, backColor=SOFT, borderPadding=10)
    para(c, "conda env create -f kemp_evb/calculation/environment.yml<br/>conda activate kemp-short<br/>python kemp_evb/calculation/scripts/run_short_evb_demo.py --check", 68, 345, 820, command_style)
    para(c, "Submit the calculation from the same repository root. Replace the example Conda path with the path on your cluster:", 68, 240, 820, BODY)
    para(c, 'CONDA_BASE="/cluster/path/to/miniforge3" sbatch kemp_evb/calculation/cluster/short_evb_demo.slurm', 68, 173, 820, command_style)
    para(c, "This bundled run evaluates a reduced five-atom path. For full-enzyme short and longer production commands, scheduler settings, and resume instructions, see the online cluster guide.", 68, 98, 820, SMALL)
    c.linkURL("https://github.com/ND7996/Kemp-EVB-Workshop/blob/main/kemp_evb/guides/CLUSTER_QUICKSTART.md", (68, 48, 850, 70), relative=0)
    c.save()
    from add_slide_notes import patch as add_slide_notes
    add_slide_notes(OUTPUT)
    from update_reference_appendix import refresh
    refresh()
    print(OUTPUT)


if __name__ == "__main__":
    build()
