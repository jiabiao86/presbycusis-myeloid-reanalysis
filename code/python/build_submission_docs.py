#!/usr/bin/env python3
"""Build the four initial-submission DOCX files from Markdown sources."""

from __future__ import annotations

import argparse
import re
from datetime import date
from pathlib import Path

from docx import Document
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


BLACK = RGBColor(0, 0, 0)
HEADER_FILL = "EDEDED"
ALT_FILL = "F7F7F7"
BORDER = "D9D9D9"
MANUSCRIPT_TITLE = (
    "Intersection of Myeloid-Associated Genes in Aging CBA/CaJ Mouse "
    "Auditory Tissues: A Secondary Reanalysis with Cross-Strain "
    "Single-Nucleus Context"
)


def set_run_font(
    run,
    size: float,
    bold: bool | None = None,
    italic: bool | None = None,
) -> None:
    run.font.name = "Arial"
    run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), "Arial")
    run.font.size = Pt(size)
    run.font.color.rgb = BLACK
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic


def add_inline_text(
    paragraph,
    text: str,
    size: float,
    bold: bool = False,
) -> None:
    """Render simple Markdown bold, italic, code, and superscript spans."""
    pattern = re.compile(r"(\*\*.+?\*\*|\*[^*]+?\*|`[^`]+`|<sup>.+?</sup>)")
    position = 0
    for match in pattern.finditer(text):
        if match.start() > position:
            run = paragraph.add_run(text[position : match.start()])
            set_run_font(run, size, bold)
        token = match.group(0)
        if token.startswith("**"):
            run = paragraph.add_run(token[2:-2])
            set_run_font(run, size, True)
        elif token.startswith("`"):
            run = paragraph.add_run(token[1:-1])
            set_run_font(run, size, bold)
        elif token.startswith("<sup>"):
            run = paragraph.add_run(token[5:-6])
            set_run_font(run, size, bold)
            run.font.superscript = True
        else:
            run = paragraph.add_run(token[1:-1])
            set_run_font(run, size, bold, True)
        position = match.end()
    if position < len(text):
        run = paragraph.add_run(text[position:])
        set_run_font(run, size, bold)


def remove_paragraph_borders(paragraph) -> None:
    paragraph_properties = paragraph._p.get_or_add_pPr()
    borders = paragraph_properties.find(qn("w:pBdr"))
    if borders is not None:
        paragraph_properties.remove(borders)


def remove_style_borders(style) -> None:
    style_properties = style.element.find(qn("w:pPr"))
    if style_properties is None:
        return
    borders = style_properties.find(qn("w:pBdr"))
    if borders is not None:
        style_properties.remove(borders)


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instruction = OxmlElement("w:instrText")
    instruction.set(qn("xml:space"), "preserve")
    instruction.text = " PAGE "
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.append(begin)
    run._r.append(instruction)
    run._r.append(end)
    set_run_font(run, 9)


def configure_document(
    doc: Document,
    document_title: str,
    body_size: float = 11,
) -> None:
    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(0.75)
    section.bottom_margin = Inches(0.75)
    section.left_margin = Inches(0.85)
    section.right_margin = Inches(0.85)
    section.header_distance = Inches(0.3)
    section.footer_distance = Inches(0.3)

    styles = doc.styles
    normal = styles["Normal"]
    normal.font.name = "Arial"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
    normal.font.size = Pt(body_size)
    normal.font.color.rgb = BLACK
    normal.paragraph_format.line_spacing = 1.15
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    title = styles["Title"]
    title.font.name = "Arial"
    title._element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
    title.font.size = Pt(20)
    title.font.bold = True
    title.font.color.rgb = BLACK
    title.paragraph_format.space_after = Pt(12)
    title.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
    remove_style_borders(title)

    for style_name, size in (
        ("Heading 1", 14),
        ("Heading 2", 12.5),
        ("Heading 3", 11.5),
    ):
        style = styles[style_name]
        style.font.name = "Arial"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = BLACK
        style.paragraph_format.keep_with_next = True
        style.paragraph_format.space_before = Pt(10)
        style.paragraph_format.space_after = Pt(5)

    add_page_number(section.footer.paragraphs[0])

    properties = doc.core_properties
    properties.title = document_title
    properties.subject = "Audited auditory transcriptome reanalysis"
    properties.author = "Jiabiao Ji and coauthors"
    properties.keywords = (
        "presbycusis; macrophage; microglia; cochlea; inferior colliculus; "
        "limma; WGCNA; CellChat"
    )


def add_title(doc: Document, text: str, size: float = 20) -> None:
    paragraph = doc.add_paragraph(style="Title")
    remove_paragraph_borders(paragraph)
    add_inline_text(paragraph, text, size, True)


def set_cell_shading(cell, fill: str) -> None:
    cell_properties = cell._tc.get_or_add_tcPr()
    shading = cell_properties.find(qn("w:shd"))
    if shading is None:
        shading = OxmlElement("w:shd")
        cell_properties.append(shading)
    shading.set(qn("w:fill"), fill)


def set_table_borders(table) -> None:
    table_properties = table._tbl.tblPr
    borders = table_properties.first_child_found_in("w:tblBorders")
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        table_properties.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        element = borders.find(qn(f"w:{edge}"))
        if element is None:
            element = OxmlElement(f"w:{edge}")
            borders.append(element)
        element.set(qn("w:val"), "single")
        element.set(qn("w:sz"), "4")
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), BORDER)


def set_cell_margins(
    cell,
    top: int = 80,
    start: int = 90,
    bottom: int = 80,
    end: int = 90,
) -> None:
    cell_properties = cell._tc.get_or_add_tcPr()
    cell_margins = cell_properties.first_child_found_in("w:tcMar")
    if cell_margins is None:
        cell_margins = OxmlElement("w:tcMar")
        cell_properties.append(cell_margins)
    for margin, value in (
        ("top", top),
        ("start", start),
        ("bottom", bottom),
        ("end", end),
    ):
        node = cell_margins.find(qn(f"w:{margin}"))
        if node is None:
            node = OxmlElement(f"w:{margin}")
            cell_margins.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def repeat_table_header(row) -> None:
    row_properties = row._tr.get_or_add_trPr()
    header = OxmlElement("w:tblHeader")
    header.set(qn("w:val"), "true")
    row_properties.append(header)


def parse_table(lines: list[str]) -> list[list[str]]:
    rows: list[list[str]] = []
    for line in lines:
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if cells and all(
            re.fullmatch(r":?-{3,}:?", cell.replace(" ", "")) for cell in cells
        ):
            continue
        rows.append(cells)
    return rows


def set_table_column_widths(table, rows: list[list[str]]) -> None:
    column_count = max(len(row) for row in rows)
    header = rows[0]
    if column_count == 3:
        widths = [1.7, 2.5, 2.6]
    elif column_count == 5 and header and header[0] == "Dataset and provenance":
        widths = [1.12, 0.96, 1.42, 1.15, 2.15]
    elif column_count == 5 and header and header[0] == "Tissue":
        widths = [1.20, 1.35, 0.55, 0.55, 3.15]
    elif column_count == 6 and header and header[0] == "Dataset":
        widths = [0.85, 1.15, 0.65, 1.45, 1.45, 1.25]
    elif column_count == 6 and header and header[0] == "Gene":
        widths = [0.65, 1.1, 1.25, 1.2, 1.25, 1.35]
    elif column_count == 7 and header and header[0] == "Tissue":
        widths = [0.95, 0.7, 0.65, 0.75, 0.75, 0.8, 2.2]
    elif column_count == 7 and header and header[0] == "Gene":
        widths = [0.72, 1.18, 1.18, 0.84, 0.84, 1.00, 1.04]
    else:
        widths = [6.8 / column_count] * column_count
    for row in table.rows:
        for cell, width in zip(row.cells, widths):
            cell.width = Inches(width)


def add_table(
    doc: Document,
    rows: list[list[str]],
    caption: str | None,
    font_size: float = 8.7,
) -> None:
    if caption:
        paragraph = doc.add_paragraph()
        paragraph.paragraph_format.keep_with_next = True
        paragraph.paragraph_format.space_before = Pt(6)
        paragraph.paragraph_format.space_after = Pt(4)
        add_inline_text(paragraph, caption, 9.3, True)

    column_count = max(len(row) for row in rows)
    table = doc.add_table(rows=0, cols=column_count)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    set_table_borders(table)

    for row_index, row in enumerate(rows):
        cells = table.add_row().cells
        for column_index in range(column_count):
            value = row[column_index] if column_index < len(row) else ""
            cell = cells[column_index]
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            set_cell_margins(cell)
            paragraph = cell.paragraphs[0]
            paragraph.paragraph_format.space_after = Pt(0)
            paragraph.paragraph_format.line_spacing = 1.0
            paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
            add_inline_text(paragraph, value, font_size, row_index == 0)
            if row_index == 0:
                set_cell_shading(cell, HEADER_FILL)
            elif row_index % 2 == 0:
                set_cell_shading(cell, ALT_FILL)
        if row_index == 0 and hasattr(table.rows[0], "_tr"):
            repeat_table_header(table.rows[0])

    set_table_column_widths(table, rows)
    spacer = doc.add_paragraph()
    spacer.paragraph_format.space_after = Pt(5)


def add_picture(doc: Document, image_path: Path, alt_text: str) -> None:
    if not image_path.exists():
        raise FileNotFoundError(image_path)
    doc.add_picture(str(image_path), width=Inches(6.55))
    paragraph = doc.paragraphs[-1]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.keep_with_next = True
    for properties in paragraph._p.xpath(".//wp:docPr"):
        properties.set("descr", alt_text)
        properties.set("title", alt_text)


def parse_markdown_table(lines: list[str], start: int) -> tuple[list[str], int]:
    table_lines: list[str] = []
    index = start
    while index < len(lines) and lines[index].strip().startswith("|"):
        table_lines.append(lines[index])
        index += 1
    return parse_table(table_lines), index


def build_manuscript(
    markdown_path: Path,
    output_path: Path,
    figure_dir: Path,
) -> None:
    lines = markdown_path.read_text(encoding="utf-8").splitlines()
    doc = Document()
    configure_document(doc, MANUSCRIPT_TITLE)

    pending_table_caption: str | None = None
    figure_caption_count = 0
    in_figures = False
    in_tables = False
    in_references = False
    index = 0

    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue

        image_match = re.fullmatch(r"!\[[^\]]*\]\(([^)]+)\)", line)
        if image_match:
            source_name = Path(image_match.group(1)).name
            image_path = figure_dir / source_name
            alt_match = re.match(r"!\[([^\]]*)\]", line)
            add_picture(
                doc,
                image_path,
                alt_match.group(1) if alt_match else source_name,
            )
            spacer = doc.add_paragraph()
            spacer.paragraph_format.space_after = Pt(4)
            index += 1
            continue

        if line.startswith("|"):
            rows, index = parse_markdown_table(lines, index)
            add_table(
                doc,
                rows,
                pending_table_caption,
                font_size=8.7,
            )
            pending_table_caption = None
            continue

        if line.startswith("# "):
            add_title(doc, line[2:].strip())
        elif line.startswith("## "):
            heading = line[3:].strip()
            in_figures = heading == "Figures"
            in_tables = heading == "Tables"
            in_references = heading == "References"
            if in_figures or in_tables:
                doc.add_page_break()
            paragraph = doc.add_paragraph(style="Heading 1")
            add_inline_text(paragraph, heading, 14, True)
        elif line.startswith("### "):
            paragraph = doc.add_paragraph(style="Heading 2")
            add_inline_text(paragraph, line[4:].strip(), 12.5, True)
        elif line.startswith("- "):
            paragraph = doc.add_paragraph(style="List Bullet")
            paragraph.paragraph_format.space_after = Pt(3)
            add_inline_text(paragraph, line[2:].strip(), 11)
        elif re.match(r"^\d+\.\s", line):
            paragraph = doc.add_paragraph()
            paragraph.paragraph_format.left_indent = Inches(0.22)
            paragraph.paragraph_format.first_line_indent = Inches(-0.22)
            paragraph.paragraph_format.space_after = Pt(4)
            add_inline_text(paragraph, line, 10.5)
        elif re.match(
            r"\*\*(?:Supplementary Figure S\d+|Figure \d+)\.",
            line,
        ):
            figure_caption_count += 1
            if in_figures and figure_caption_count > 1:
                doc.add_page_break()
            paragraph = doc.add_paragraph()
            paragraph.paragraph_format.keep_with_next = True
            paragraph.paragraph_format.space_before = Pt(4)
            paragraph.paragraph_format.space_after = Pt(5)
            add_inline_text(paragraph, line, 10.2, True)
        elif line.startswith("**Table ") and line.endswith("**"):
            if line.startswith(("**Table 3.", "**Table 5.")):
                doc.add_page_break()
            pending_table_caption = line.strip("*")
        else:
            paragraph = doc.add_paragraph()
            if in_references:
                paragraph.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
                paragraph.paragraph_format.line_spacing = 1.0
                paragraph.paragraph_format.space_after = Pt(4)
                add_inline_text(paragraph, line, 10)
            else:
                add_inline_text(paragraph, line, 11)
            if line.startswith("**Keywords:**"):
                paragraph.add_run().add_break(WD_BREAK.PAGE)
        index += 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(output_path)


def markdown_blocks(markdown_path: Path) -> list[str]:
    lines = markdown_path.read_text(encoding="utf-8").splitlines()
    blocks: list[str] = []
    current: list[str] = []
    for line in lines:
        if line.strip():
            current.append(line.rstrip())
        elif current:
            blocks.append("\n".join(current))
            current = []
    if current:
        blocks.append("\n".join(current))
    return blocks


def build_cover_letter(markdown_path: Path, output_path: Path) -> None:
    doc = Document()
    configure_document(doc, "Cover Letter", body_size=10)
    section = doc.sections[0]
    section.top_margin = Inches(0.55)
    section.bottom_margin = Inches(0.55)
    section.left_margin = Inches(0.70)
    section.right_margin = Inches(0.70)
    normal = doc.styles["Normal"]
    normal.font.size = Pt(10)
    normal.paragraph_format.line_spacing = 0.93
    normal.paragraph_format.space_after = Pt(3)
    title = doc.styles["Title"]
    title.font.size = Pt(15)
    title.paragraph_format.space_after = Pt(6)
    blocks = markdown_blocks(markdown_path)
    title_written = False
    date_written = False
    signature_block_started = False
    for block in blocks:
        stripped = block.strip()
        if stripped.startswith("# "):
            add_title(doc, stripped[2:].strip(), 15)
            title_written = True
            continue
        if title_written and not date_written:
            paragraph = doc.add_paragraph(date.today().strftime("%B %-d, %Y"))
            paragraph.paragraph_format.space_after = Pt(5)
            date_written = True
        for line in stripped.splitlines():
            paragraph = doc.add_paragraph()
            if signature_block_started:
                paragraph.paragraph_format.space_after = Pt(1)
            else:
                paragraph.paragraph_format.space_after = (
                    Pt(3) if "\n" not in stripped else Pt(0)
                )
            add_inline_text(paragraph, line, 10)
            if line.strip() == "Sincerely,":
                signature_block_started = True
    output_path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(output_path)


def build_highlights(markdown_path: Path, output_path: Path) -> None:
    doc = Document()
    configure_document(doc, "Highlights", body_size=11)
    add_title(doc, "Highlights", 18)
    for line in markdown_path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("# "):
            continue
        if stripped.startswith("- "):
            paragraph = doc.add_paragraph(style="List Bullet")
            paragraph.paragraph_format.space_after = Pt(8)
            add_inline_text(paragraph, stripped[2:].strip(), 11)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(output_path)


def build_figure_legends(
    markdown_path: Path,
    output_path: Path,
    figure_dir: Path,
) -> None:
    doc = Document()
    configure_document(doc, "Figure Legends", body_size=10.5)
    add_title(doc, "Figure Legends", 18)
    lines = markdown_path.read_text(encoding="utf-8").splitlines()
    legend_blocks: list[str] = []
    for line in lines:
        stripped = line.strip()
        if re.match(
            r"\*\*(?:Graphical abstract|Figure \d+|Supplementary Figure S\d+)\.",
            stripped,
        ):
            legend_blocks.append(stripped)

    figure_paths = {
        "Graphical abstract": "graphical_abstract_v14.png",
        **{
            f"Figure {number}": f"figure{number}_{slug}.png"
            for number, slug in {
                1: "study_audit",
                2: "frozen_limma",
                3: "intersection_sensitivity",
                4: "cellular_context",
                5: "model_identifiability",
            }.items()
        },
        **{
            f"Supplementary Figure S{number}": f"figureS{number}_{slug}.png"
            for number, slug in {
                1: "wgcna",
                2: "external_context",
                3: "cellchat_age_overview",
                4: "cellchat_focus",
                5: "myeloid_subclustering",
                6: "full_atlas",
            }.items()
        },
    }

    for block_index, block in enumerate(legend_blocks):
        if block_index:
            doc.add_page_break()
        paragraph = doc.add_paragraph()
        paragraph.paragraph_format.keep_with_next = True
        paragraph.paragraph_format.space_after = Pt(8)
        add_inline_text(paragraph, block, 10.5)
        label_match = re.match(
            r"\*\*(Graphical abstract|Figure \d+|Supplementary Figure S\d+)\.",
            block,
        )
        if label_match:
            image_path = figure_dir / figure_paths[label_match.group(1)]
            add_picture(doc, image_path, label_match.group(1))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(output_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        help="Folder containing the manuscript and supporting Markdown sources.",
    )
    parser.add_argument(
        "--package-root",
        required=True,
        type=Path,
        help="Target submission-package root.",
    )
    parser.add_argument(
        "--figure-dir",
        required=True,
        type=Path,
        help="Folder containing the approved PNG files.",
    )
    args = parser.parse_args()

    project_root = (
        args.project_root.resolve()
        if args.project_root
        else Path(__file__).resolve().parents[3]
    )
    package_root = args.package_root.resolve()
    figure_dir = args.figure_dir.resolve()

    def source_file(relative: str) -> Path:
        direct = project_root / relative
        if direct.is_file():
            return direct
        packaged = package_root / relative
        if packaged.is_file():
            return packaged
        matches = sorted(project_root.rglob(Path(relative).name))
        if len(matches) == 1:
            return matches[0]
        raise FileNotFoundError(
            f"Could not locate {Path(relative).name} under {project_root} "
            f"or {packaged.parent}"
        )

    build_manuscript(
        source_file("Manuscript/Manuscript.md"),
        package_root / "Manuscript" / "Peripheral_Central_Overlap_Manuscript.docx",
        figure_dir,
    )
    build_cover_letter(
        source_file("Cover_Letter/Cover_Letter.md"),
        package_root / "Cover_Letter" / "Cover_Letter.docx",
    )
    build_highlights(
        source_file("Highlights/Highlights.md"),
        package_root / "Highlights" / "Highlights.docx",
    )
    build_figure_legends(
        source_file("Figure_Legends/Figure_Legends.md"),
        package_root / "Figure_Legends" / "Figure_Legends.docx",
        figure_dir,
    )

    print(f"Built four initial-submission DOCX files under {package_root}")


if __name__ == "__main__":
    main()
