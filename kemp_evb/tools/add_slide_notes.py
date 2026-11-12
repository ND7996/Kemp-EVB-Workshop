"""Finish the slide PDF with the custom mechanism and labeled profile."""
from pathlib import Path

import pymupdf as fitz


def patch(pdf_path: Path) -> None:
    pdf = fitz.open(pdf_path)

    # Remove the welcome paragraph as requested.
    page = pdf[0]
    page.add_redact_annot(fitz.Rect(68, 145, 850, 280), fill=(1, 1, 1))
    page.apply_redactions()

    # Refresh the scheme image after correcting the product C≡N bond.
    page = pdf[1]
    for image in page.get_images(full=True):
        xref = image[0]
        if page.get_image_rects(xref):
            page.replace_image(xref, filename=str(pdf_path.parent / "assets" / "kemp_states.png"))
            break

    # Keep the atom key/output panel, replacing it first so this is repeatable.
    page = pdf[2]
    page.add_redact_annot(fitz.Rect(68, 428, 892, 490), fill=(0.918, 0.953, 0.949))
    page.apply_redactions()
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

    # Replace the explanatory cards with the profile itself; its axes carry
    # the labels and units.
    page = pdf[7]
    page.add_redact_annot(fitz.Rect(60, 88, 900, 475), fill=(1, 1, 1))
    page.apply_redactions()
    page.insert_text((68, 128), "The saved reaction profile", fontsize=28, fontname="hebo", color=(0.12, 0.18, 0.21))
    page.insert_image(
        fitz.Rect(175, 140, 785, 470),
        filename=str(pdf_path.parent / "assets" / "source_profile_slide.png"),
        keep_proportion=True,
    )
    page.insert_text((68, 488), "Saved short-run profile; not validated.", fontsize=10, fontname="helv", color=(0.32, 0.39, 0.42))

    # Remove the obsolete explanatory paragraph from the closing slide.
    page = pdf[9]
    page.add_redact_annot(fitz.Rect(68, 345, 890, 380), fill=(1, 1, 1))
    page.apply_redactions()

    patched = Path(pdf_path).with_name(Path(pdf_path).stem + ".patched.pdf")
    pdf.save(patched, garbage=4, deflate=True)
    pdf.close()
    patched.replace(pdf_path)


if __name__ == "__main__":
    patch(Path(__file__).resolve().parents[1] / "Kemp_EVB_workshop.pdf")
