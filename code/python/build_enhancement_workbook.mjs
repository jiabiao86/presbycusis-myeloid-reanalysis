import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis/manuscript/supplementary_package";
const manifest = JSON.parse(await fs.readFile(path.join(root, "enhancement_tables_data.json"), "utf8"));
const part = process.argv[2] ?? "cellchat";
const partConfig = {
  cellchat: {
    prefixes: ["S19", "S20"],
    fileName: "Supplementary_Tables_S19-S20_CellChat.xlsx",
  },
  enrichment: {
    prefixes: ["S21", "S22", "S23", "S24"],
    fileName: "Supplementary_Tables_S21-S24_Enrichment.xlsx",
  },
};
if (!partConfig[part]) {
  throw new Error(`Unknown enhancement workbook part: ${part}`);
}
const workbook = Workbook.create();
const contents = workbook.worksheets.add("Contents");
contents.showGridLines = false;
contents.tabColor = "#1F4E79";

function sheetSpecs() {
  const specs = [];
  for (const table of manifest.tables) {
    if (table.key === "S19_CellChatDB_communication_scores") {
      [
        [3, "S19A_communication_3M"],
        [12, "S19B_communication_12M"],
        [24, "S19C_communication_24M"],
      ].forEach(([age, sheetName]) => {
        const rows = table.rows.filter((row) => Number(row.age_months) === Number(age));
        specs.push({
          tableNumber: `S19${sheetName.split("_")[0].slice(-1)}`,
          title: `${table.title} ${age} months`,
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
      title: table.title,
      columns: table.columns,
      rows: table.rows,
      row_count: table.row_count,
      column_count: table.column_count,
      sheetName: table.key.slice(0, 31),
    });
  }
  return specs;
}

function columnWidth(columnName) {
  const name = columnName.toLowerCase();
  if (name === "gene") return 18;
  if (name === "term") return 58;
  if (name.includes("pathway")) return 22;
  if (name.includes("sender") || name.includes("receiver") || name.includes("cell_type")) return 24;
  if (name.includes("interaction") || name.includes("complex")) return 28;
  if (name.includes("p_value") || name.includes("fdr")) return 16;
  if (name.includes("score") || name.includes("mean") || name.includes("percent")) return 16;
  return Math.min(24, Math.max(11, columnName.length + 2));
}

function formatForColumn(columnName) {
  const name = columnName.toLowerCase();
  if (/p_value|fdr/.test(name)) return "0.00E+00";
  if (/score|mean|percent|z_score/.test(name)) return "0.000";
  if (/age_months|rank/.test(name)) return "0";
  return "@";
}

function styleHeader(sheet, columnCount) {
  const range = sheet.getRangeByIndexes(0, 0, 1, columnCount);
  range.format.fill = "#1F4E79";
  range.format.font = { bold: true, color: "#FFFFFF", size: 10, typeface: "Arial" };
  range.format.horizontalAlignment = "center";
  range.format.verticalAlignment = "center";
  range.format.wrapText = true;
  range.format.rowHeight = 26;
  range.format.borders = { preset: "all", style: "thin", color: "#D9D9D9" };
}

const allSpecs = sheetSpecs();
const specs = allSpecs.filter((spec) =>
  partConfig[part].prefixes.includes(spec.tableNumber.slice(0, 3))
);
contents.getRange("A1").values = [[manifest.title]];
contents.getRange("A2").values = [[manifest.manuscript_title]];
contents.getRange("A4:E4").values = [["Table", "Title", "Rows", "Columns", "Worksheet"]];
contents.getRangeByIndexes(4, 0, specs.length, 5).values = specs.map((spec) => [
  spec.tableNumber,
  spec.title,
  spec.row_count,
  spec.column_count,
  spec.sheetName,
]);
contents.getRange("A1").format.font = { bold: true, size: 16, typeface: "Arial" };
contents.getRange("A2").format.font = { italic: true, size: 10, typeface: "Arial" };
contents.getRange("A4:E4").format.fill = "#1F4E79";
contents.getRange("A4:E4").format.font = { bold: true, color: "#FFFFFF", typeface: "Arial" };
contents.getRangeByIndexes(4, 0, specs.length, 5).format.borders = {
  preset: "all",
  style: "thin",
  color: "#D9D9D9",
};
contents.getRangeByIndexes(4, 0, specs.length, 5).format.font = { typeface: "Arial", size: 10 };
contents.getRange("A:A").format.columnWidth = 12;
contents.getRange("B:B").format.columnWidth = 76;
contents.getRange("C:D").format.columnWidth = 12;
contents.getRange("E:E").format.columnWidth = 28;
contents.freezePanes.freezeRows(4);

for (const spec of specs) {
  const sheet = workbook.worksheets.add(spec.sheetName);
  sheet.showGridLines = false;
  sheet.tabColor = "#8C6BB1";
  const matrix = [
    spec.columns,
    ...spec.rows.map((row) => spec.columns.map((column) => row[column] ?? null)),
  ];
  sheet.getRangeByIndexes(0, 0, matrix.length, spec.column_count).values = matrix;
  styleHeader(sheet, spec.column_count);
  for (let columnIndex = 0; columnIndex < spec.columns.length; columnIndex += 1) {
    const columnName = spec.columns[columnIndex];
    if (spec.row_count <= 5000) {
      sheet
        .getRangeByIndexes(1, columnIndex, spec.row_count, 1)
        .setNumberFormat(formatForColumn(columnName));
    }
    sheet.getRangeByIndexes(0, columnIndex, 1, 1).format.columnWidth = columnWidth(columnName);
  }
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(1);
}

workbook.recalculate();
const result = await workbook.inspect({
  kind: "table",
  sheetId: "Contents",
  range: `A1:E${Math.max(6, specs.length + 5)}`,
  include: "values,formulas",
  tableMaxRows: 20,
  tableMaxCols: 8,
});
console.log(result.ndjson);
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 100 },
  summary: "enhancement workbook formula error scan",
});
console.log(errors.ndjson);
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(path.join(root, partConfig[part].fileName));
console.log(`Wrote ${partConfig[part].fileName}`);
