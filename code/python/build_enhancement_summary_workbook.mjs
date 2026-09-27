import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis/manuscript/supplementary_package";
const manifest = JSON.parse(await fs.readFile(path.join(root, "enhancement_tables_data.json"), "utf8"));
const workbook = Workbook.create();
const contents = workbook.worksheets.add("Contents");
contents.showGridLines = false;
contents.tabColor = "#1F4E79";

const previewTables = manifest.tables.map((table) => ({
  ...table,
  previewRows: table.rows.slice(0, 100),
}));

contents.getRange("A1").values = [["Supplementary Tables S19-S25"]];
contents.getRange("A2").values = [[manifest.manuscript_title]];
contents.getRange("A4:E4").values = [["Table", "Title", "Complete rows", "Columns", "CSV data file"]];
contents.getRangeByIndexes(4, 0, previewTables.length, 5).values = previewTables.map((table) => [
  table.key.split("_")[0],
  table.title,
  table.row_count,
  table.column_count,
  `${table.key.replaceAll("_", "-")}.csv`,
]);
contents.getRange("A1").format.font = { bold: true, size: 16, typeface: "Arial" };
contents.getRange("A2").format.font = { italic: true, size: 10, typeface: "Arial" };
contents.getRange("A4:E4").format.fill = "#1F4E79";
contents.getRange("A4:E4").format.font = { bold: true, color: "#FFFFFF", size: 10, typeface: "Arial" };
contents.getRangeByIndexes(4, 0, previewTables.length, 5).format.borders = {
  preset: "all",
  style: "thin",
  color: "#D9D9D9",
};
contents.getRangeByIndexes(4, 0, previewTables.length, 5).format.font = { typeface: "Arial", size: 10 };
contents.getRange("A:A").format.columnWidth = 12;
contents.getRange("B:B").format.columnWidth = 72;
contents.getRange("C:D").format.columnWidth = 14;
contents.getRange("E:E").format.columnWidth = 44;
contents.freezePanes.freezeRows(4);

for (const table of previewTables) {
  const sheet = workbook.worksheets.add(table.key.slice(0, 31));
  sheet.showGridLines = false;
  sheet.tabColor = "#A9C4E2";
  const matrix = [
    table.columns,
    ...table.previewRows.map((row) => table.columns.map((column) => row[column] ?? null)),
  ];
  sheet.getRangeByIndexes(0, 0, matrix.length, table.column_count).values = matrix;
  const header = sheet.getRangeByIndexes(0, 0, 1, table.column_count);
  header.format.fill = "#1F4E79";
  header.format.font = { bold: true, color: "#FFFFFF", size: 10, typeface: "Arial" };
  header.format.horizontalAlignment = "center";
  header.format.verticalAlignment = "center";
  header.format.wrapText = true;
  header.format.rowHeight = 26;
  header.format.borders = { preset: "all", style: "thin", color: "#D9D9D9" };
  for (let index = 0; index < table.columns.length; index += 1) {
    const name = table.columns[index];
    sheet.getRangeByIndexes(0, index, 1, 1).format.columnWidth =
      name === "term" ? 58 : name.includes("gene") ? 20 : Math.min(28, Math.max(12, name.length + 2));
  }
  sheet.freezePanes.freezeRows(1);
}

workbook.recalculate();
const preview = await workbook.render({ sheetName: "Contents", scale: 1, format: "png", autoCrop: "all" });
await fs.writeFile(
  path.join(root, "Supplementary_Tables_S19-S24_Summary_Contents.png"),
  new Uint8Array(await preview.arrayBuffer())
);
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(path.join(root, "Supplementary_Tables_S19-S25_Summary.xlsx"));
console.log("Wrote Supplementary_Tables_S19-S25_Summary.xlsx");
