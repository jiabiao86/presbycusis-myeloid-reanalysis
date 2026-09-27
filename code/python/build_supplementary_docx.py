#!/usr/bin/env python3
"""Build the human-readable supplementary table index and captions."""

from __future__ import annotations

import json
from pathlib import Path

from docx import Document
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis/manuscript/supplementary_package")
MANIFESTS = [
    ROOT / "supplementary_tables_data.json",
    ROOT / "enhancement_tables_data.json",
    ROOT / "posthoc_tables_data.json",
    ROOT / "full_cellchat_tables_data.json",
]
OUTPUT = ROOT / "Supplementary_Tables_S1-S37.docx"

BLACK = RGBColor(0, 0, 0)
WHITE = RGBColor(255, 255, 255)
HEADER_FILL = "1F4E79"
ALT_FILL = "EAF2F8"
BORDER = "D9D9D9"


def set_run_font(run, size: float, bold: bool | None = None) -> None:
    run.font.name = "Arial"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
    run.font.size = Pt(size)
    run.font.color.rgb = BLACK
    if bold is not None:
        run.bold = bold


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_table_borders(table) -> None:
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.first_child_found_in("w:tblBorders")
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        element = borders.find(qn(f"w:{edge}"))
        if element is None:
            element = OxmlElement(f"w:{edge}")
            borders.append(element)
        element.set(qn("w:val"), "single")
        element.set(qn("w:sz"), "4")
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), BORDER)


def set_cell_margins(cell, top=90, start=110, bottom=90, end=110) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for margin, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{margin}"))
        if node is None:
            node = OxmlElement(f"w:{margin}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


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


def configure_document(doc: Document) -> None:
    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(0.75)
    section.bottom_margin = Inches(0.75)
    section.left_margin = Inches(0.85)
    section.right_margin = Inches(0.85)

    normal = doc.styles["Normal"]
    normal.font.name = "Arial"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
    normal.font.size = Pt(10.5)
    normal.font.color.rgb = BLACK
    normal.paragraph_format.line_spacing = 1.12
    normal.paragraph_format.space_after = Pt(5)
    normal.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    title = doc.styles["Title"]
    title.font.name = "Arial"
    title._element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
    title.font.size = Pt(20)
    title.font.bold = True
    title.font.color.rgb = BLACK
    title.paragraph_format.space_after = Pt(14)

    for style_name, size in (("Heading 1", 14), ("Heading 2", 12), ("Heading 3", 11)):
        style = doc.styles[style_name]
        style.font.name = "Arial"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = BLACK
        style.paragraph_format.keep_with_next = True
        style.paragraph_format.space_before = Pt(10)
        style.paragraph_format.space_after = Pt(4)

    add_page_number(section.footer.paragraphs[0])
    props = doc.core_properties
    props.title = "Supplementary Tables"
    props.subject = "Peripheral-Central Convergence of Macrophage and Microglial Programs in Presbycusis Severity"
    props.author = "Manuscript supplement"
    props.keywords = "supplementary tables; presbycusis; transcriptomics; WGCNA"


def add_table(doc: Document, columns: list[str], rows: list[dict], preview_rows: int | None = None) -> None:
    data_rows = rows if preview_rows is None else rows[:preview_rows]
    table = doc.add_table(rows=0, cols=len(columns))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    set_table_borders(table)

    for row_index, values in enumerate([columns] + [[row.get(column) for column in columns] for row in data_rows]):
        cells = table.add_row().cells
        for column_index, value in enumerate(values):
            cell = cells[column_index]
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            set_cell_margins(cell)
            paragraph = cell.paragraphs[0]
            paragraph.paragraph_format.space_after = Pt(0)
            paragraph.paragraph_format.line_spacing = 1.0
            text = "" if value is None else str(value)
            run = paragraph.add_run(text)
            set_run_font(run, 8.2, row_index == 0)
            if row_index == 0:
                set_cell_shading(cell, HEADER_FILL)
                run.font.color.rgb = WHITE
            elif row_index % 2 == 0:
                set_cell_shading(cell, ALT_FILL)
        if row_index == 0:
            repeat_table_header(table.rows[0])
    spacer = doc.add_paragraph()
    spacer.paragraph_format.space_after = Pt(4)


def data_source_for_key(key: str) -> str:
    number = int(key.split("_")[0].lstrip("S"))
    if number <= 11:
        return "Supplementary_Tables_S1-S11.xlsx"
    if number <= 18:
        return "Supplementary_Tables_S12-S18.xlsx"
    filename = key.replace("_", "-") + ".csv"
    if number <= 25:
        return f"enhancement_tables/{filename}"
    if number <= 34:
        return f"posthoc_tables/{filename}"
    return f"full_cellchat_tables/{filename}"


def worksheet_for_key(key: str) -> str:
    mapping = {
        "S7_WGCNA_module_traits_preservation": "S7_WGCNA_modules",
        "S9_GSE274279_candidate_expression": "S9_candidate_expression",
        "S11_Peripheral_central_comparison": "S11_Peripheral_central",
        "S13_GSE153882_all_genes": "S13A_GSE153882_IHC; S13B_GSE153882_OHC",
    }
    if key == "S3_Limma_all_results":
        return "S3A_MA_vs_YC; S3B_MP_vs_MA; S3C_MP_vs_YC; S3D_SP_vs_MA; S3E_SP_vs_MP; S3F_SP_vs_YC"
    if key.startswith(("S19_", "S20_", "S21_", "S22_", "S23_", "S24_")):
        return "CSV data file"
    if key.startswith(("S35_", "S36_", "S37_")):
        return "CSV data file"
    return mapping.get(key, key[:31])


def display_columns(columns: list[str]) -> str:
    return ", ".join(columns)


def main() -> None:
    manifests = [json.loads(path.read_text(encoding="utf-8")) for path in MANIFESTS]
    tables = [table for manifest in manifests for table in manifest["tables"]]
    doc = Document()
    configure_document(doc)

    title = doc.add_paragraph(style="Title")
    run = title.add_run("Supplementary Tables")
    set_run_font(run, 20, True)

    subtitle = doc.add_paragraph()
    run = subtitle.add_run(manifests[0]["manuscript_title"])
    set_run_font(run, 11, True)
    subtitle.paragraph_format.space_after = Pt(10)

    overview = doc.add_paragraph()
    overview.add_run(
        "This supplement contains caption pages and structured data indexes for Supplementary Tables S1-S37. "
        "Large machine-readable result tables are provided in the accompanying Excel workbooks. "
        "Each table section states the workbook, worksheet, row count, and column names."
    )

    doc.add_heading("Table Index", level=1)
    index_rows = []
    for table in tables:
        key = table["key"]
        index_rows.append(
            {
                "Table": key.split("_")[0],
                "Title": table["title"],
                "Rows": table["row_count"],
                "Columns": table["column_count"],
                "Data file": data_source_for_key(key),
                "Worksheet": worksheet_for_key(key),
            }
        )
    add_table(doc, ["Table", "Title", "Rows", "Columns", "Data file", "Worksheet"], index_rows)

    for table in tables:
        doc.add_page_break()
        heading = doc.add_heading(table["title"], level=1)
        for run in heading.runs:
            set_run_font(run, 14, True)

        description = doc.add_paragraph()
        description.add_run(table["description"])

        source = doc.add_paragraph()
        run = source.add_run(
            f"Data file: {data_source_for_key(table['key'])}. "
            f"Worksheet(s): {worksheet_for_key(table['key'])}. "
            f"Rows: {table['row_count']:,}. Columns: {table['column_count']}."
        )
        set_run_font(run, 9.5)

        if table.get("notes"):
            note = doc.add_paragraph()
            run = note.add_run(f"Note: {table['notes']}")
            set_run_font(run, 9.5)

        columns = doc.add_paragraph()
        run = columns.add_run(f"Columns: {display_columns(table['columns'])}")
        set_run_font(run, 9.5)

        if table["row_count"] <= 25 and table["column_count"] <= 10:
            add_table(doc, table["columns"], table["rows"])
        elif table["row_count"] <= 600 and table["column_count"] <= 10:
            add_table(doc, table["columns"], table["rows"], preview_rows=12)
            preview_note = doc.add_paragraph()
            run = preview_note.add_run(
                "The complete table is provided in the accompanying Excel workbook."
            )
            set_run_font(run, 9.5)
        else:
            preview_note = doc.add_paragraph()
            run = preview_note.add_run(
                "This table is too large for readable presentation in Word. "
                "The complete machine-readable table is provided in the accompanying Excel workbook."
            )
            set_run_font(run, 9.5)

    doc.save(OUTPUT)
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
