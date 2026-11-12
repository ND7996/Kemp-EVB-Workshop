"""Add the atom-label key and method outputs to the slide PDF."""
from pathlib import Path

import pymupdf as fitz


def patch(pdf_path: Path) -> None:
    pdf = fitz.open(pdf_path)
    page = pdf[2]
    page.draw_rect(fitz.Rect(68, 428, 892, 490), color=None, fill=(0.918, 0.953, 0.949))
    page.insert_text((82, 445), "ATOM LABELS", fontsize=10, fontname="hebo", color=(0.09, 0.42, 0.45))
    page.insert_textbox(
        fitz.Rect(82, 451, 475, 484),
        "C3 is the substrate carbon that loses Ht; Ht is the transferred hydrogen; O1 and N2 are the ring oxygen and nitrogen. Asp162 is enzyme aspartate 162; OD1 is its accepting side-chain oxygen.",
        fontsize=9, fontname="helv", color=(0.12, 0.18, 0.21), lineheight=1.15,
    )
    page.insert_text((495, 445), "WHAT THE METHODS GIVE", fontsize=10, fontname="hebo", color=(0.09, 0.42, 0.45))
    page.insert_textbox(
        fitz.Rect(495, 451, 878, 488),
        "Q6 reports activation-barrier estimates for enzyme variants. A converged OpenMM run gives a free-energy profile, barrier, and reaction free energy. The saved profile is unvalidated.",
        fontsize=8.5, fontname="helv", color=(0.12, 0.18, 0.21), lineheight=1.12,
    )
    patched = Path(pdf_path).with_name(Path(pdf_path).stem + ".patched.pdf")
    pdf.save(patched, garbage=4, deflate=True)
    pdf.close()
    patched.replace(pdf_path)


if __name__ == "__main__":
    patch(Path(__file__).resolve().parents[1] / "Kemp_EVB_workshop.pdf")
