import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis/manuscript/supplementary_package";
const manifestPath = path.join(root, "supplementary_tables_data.json");
const part = process.argv[2] ?? "part1";
const manifest = JSON.parse(await fs.readFile(manifestPath, "utf8"));

const partDefinitions = {
  part1: {
    file: "Supplementary_Tables_S1-S11.xlsx",
    prefixes: ["S1_", "S2_", "S3_", "S4_", "S5_", "S6_", "S7_", "S8_", "S9_", "S10_", "S11_"],
  },
  part2: {
    file: "Supplementary_Tables_S12-S18.xlsx",
    prefixes: ["S12_", "S13_", "S14_", "S15_", "S16_", "S17_", "S18_"],
  },
};

if (!partDefinitions[part]) {
  throw new Error(`Unknown workbook part: ${part}`);
}

function sheetNameFor(key) {
  const mapping = {
    S7_WGCNA_module_traits_preservation: "S7_WGCNA_modules",
    S9_GSE274279_candidate_expression: "S9_candidate_expression",
    S11_Peripheral_central_comparison: "S11_Peripheral_central",
  };
  return mapping[key] ?? key.slice(0, 31);
}

function isNumeric(value) {
  return typeof value === "number" && Number.isFinite(value);
}

function formatForColumn(columnName, rows) {
  const name = columnName.toLowerCase();
  if (/(^|_)p_value$|^p_value$|_p$|fdr|adj\.p|zsummary_p|observed_p/.test(name)) {
    return "0.00E+00";
  }
  if (/correlation|log2fc|severity_effect|fold_change|welch_t|moderated_t|average_expression|log_odds|mean_log_expression|percent_expressing/.test(name)) {
    return "0.000";
  }
  if (/score|rank|genes_used|severity_score|age_months|cells|selected_count|row_count/.test(name)) {
    return "0";
  }
  const numeric = rows.find((row) => isNumeric(row[columnName]));
  return numeric ? "General" : "@";
}

function columnWidth(columnName) {
  const name = columnName.toLowerCase();
  if (name === "gene") return 18;
  if (name === "comparison" || name === "dataset" || name === "cell_type") return 20;
  if (name.includes("sample") || name.includes("gsm") || name.includes("accession")) return 18;
  if (name.includes("title") || name.includes("description") || name.includes("original_title")) return 42;
  if (name.includes("p_value") || name.includes("fdr") || name.includes("zsummary")) return 15;
  if (name.includes("correlation") || name.includes("log2fc") || name.includes("fold_change")) return 16;
  if (name.includes("task") || name.includes("model") || name.includes("module")) return 22;
  return Math.min(24, Math.max(11, columnName.length + 2));
}

function styleHeader(sheet, columnCount) {
  const headerRange = sheet.getRangeByIndexes(0, 0, 1, columnCount);
  headerRange.format.fill = "#1F4E79";
  headerRange.format.font = { bold: true, color: "#FFFFFF", size: 10, typeface: "Arial" };
  headerRange.format.horizontalAlignment = "center";
  headerRange.format.verticalAlignment = "center";
  headerRange.format.wrapText = true;
  headerRange.format.rowHeight = 26;
  headerRange.format.borders = { preset: "all", style: "thin", color: "#D9D9D9" };
}

function buildSheetSpecs(tables) {
  const specs = [];
  for (const table of tables) {
    if (table.key === "S3_Limma_all_results") {
      const comparisons = [...new Set(table.rows.map((row) => row.comparison))];
      comparisons.forEach((comparison, index) => {
        const suffix = String.fromCharCode(65 + index);
        const rows = table.rows.filter((row) => row.comparison === comparison);
        specs.push({
          tableNumber: `S3${suffix}`,
          key: `${table.key}_${suffix}`,
          title: `${table.title} ${comparison}`,
          columns: table.columns,
          rows,
          row_count: rows.length,
          column_count: table.column_count,
          sheetName: `S3${suffix}_${comparison}`.slice(0, 31),
        });
      });
      continue;
    }

    if (table.key === "S13_GSE153882_all_genes") {
      const groups = [
        ["A", "Inner hair cells", "S13A_GSE153882_IHC"],
        ["B", "Outer hair cells", "S13B_GSE153882_OHC"],
      ];
      groups.forEach(([suffix, cellType, sheetName]) => {
        const rows = table.rows.filter((row) => row.cell_type === cellType);
        specs.push({
          tableNumber: `S13${suffix}`,
          key: `${table.key}_${suffix}`,
          title: `${table.title} ${cellType}`,
          columns: table.columns,
          rows,
          row_count: rows.length,
          column_count: table.column_count,
          sheetName,
        });
      });
      continue;
    }

    specs.push({
      tableNumber: table.key.split("_")[0],
      key: table.key,
      title: table.title,
      columns: table.columns,
      rows: table.rows,
      row_count: table.row_count,
      column_count: table.column_count,
      sheetName: sheetNameFor(table.key),
    });
  }
  return specs;
}

function tableInPart(table) {
  return partDefinitions[part].prefixes.some((prefix) => table.key.startsWith(prefix));
}

const selectedTables = manifest.tables.filter(tableInPart);
const sheetSpecs = buildSheetSpecs(selectedTables);

const workbook = Workbook.create();
const contents = workbook.worksheets.add("Contents");
contents.showGridLines = false;
contents.tabColor = "#1F4E79";

contents.getRange("A1").values = [[manifest.title]];
contents.getRange("A2").values = [[manifest.manuscript_title]];
contents.getRange("A4:E4").values = [["Table", "Title", "Rows", "Columns", "Worksheet"]];
contents.getRangeByIndexes(4, 0, sheetSpecs.length, 5).values = sheetSpecs.map((table) => [
  table.tableNumber,
  table.title,
  table.row_count,
  table.column_count,
  table.sheetName,
]);

contents.getRange("A1").format.font = { bold: true, size: 16, color: "#000000", typeface: "Arial" };
contents.getRange("A2").format.font = { italic: true, size: 10, color: "#444444", typeface: "Arial" };
contents.getRange("A4:E4").format.fill = "#1F4E79";
contents.getRange("A4:E4").format.font = { bold: true, color: "#FFFFFF", size: 10, typeface: "Arial" };
contents.getRange("A4:E4").format.horizontalAlignment = "center";
contents.getRange("A4:E4").format.verticalAlignment = "center";
contents.getRange("A4:E4").format.rowHeight = 24;
contents.getRangeByIndexes(4, 0, sheetSpecs.length, 5).format.borders = {
  preset: "all",
  style: "thin",
  color: "#D9D9D9",
};
contents.getRangeByIndexes(4, 0, sheetSpecs.length, 5).format.verticalAlignment = "center";
contents.getRangeByIndexes(4, 0, sheetSpecs.length, 5).format.font = { typeface: "Arial", size: 10 };
contents.getRange("A:A").format.columnWidth = 12;
contents.getRange("B:B").format.columnWidth = 78;
contents.getRange("C:D").format.columnWidth = 12;
contents.getRange("E:E").format.columnWidth = 28;
contents.getRange("C:D").format.horizontalAlignment = "center";
contents.freezePanes.freezeRows(4);

for (const table of sheetSpecs) {
  const sheet = workbook.worksheets.add(table.sheetName);
  sheet.showGridLines = false;
  sheet.tabColor = table.tableNumber.startsWith("S3") ||
    table.tableNumber.startsWith("S4") ||
    table.tableNumber.startsWith("S5") ||
    table.tableNumber.startsWith("S12") ||
    table.tableNumber.startsWith("S13")
    ? "#5B9BD5"
    : "#A9C4E2";

  const matrix = [
    table.columns,
    ...table.rows.map((row) => table.columns.map((column) => row[column] ?? null)),
  ];
  sheet.getRangeByIndexes(0, 0, matrix.length, table.column_count).values = matrix;
  styleHeader(sheet, table.column_count);

  for (let columnIndex = 0; columnIndex < table.columns.length; columnIndex += 1) {
    const columnName = table.columns[columnIndex];
    if (table.row_count <= 5000) {
      const columnRange = sheet.getRangeByIndexes(1, columnIndex, table.row_count, 1);
      columnRange.setNumberFormat(formatForColumn(columnName, table.rows));
    }
    sheet.getRangeByIndexes(0, columnIndex, 1, 1).format.columnWidth = columnWidth(columnName);
  }

  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(1);
}

workbook.recalculate();

const contentsCheck = await workbook.inspect({
  kind: "table",
  sheetId: "Contents",
  range: `A1:E${Math.max(6, sheetSpecs.length + 5)}`,
  include: "values,formulas",
  tableMaxRows: 25,
  tableMaxCols: 8,
});
console.log(contentsCheck.ndjson);

const formulaErrors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(formulaErrors.ndjson);

const outputPath = path.join(root, partDefinitions[part].file);
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(`Wrote ${outputPath}`);
