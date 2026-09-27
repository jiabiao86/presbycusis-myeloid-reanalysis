import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const root = "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis/manuscript/supplementary_package";
const files = [
  ["Supplementary_Tables_S1-S11.xlsx", ["Contents", "S1_Sample_metadata", "S3A_MA_vs_YC", "S11_Peripheral_central"]],
  ["Supplementary_Tables_S12-S18.xlsx", ["Contents", "S12_GSE233798_all_genes", "S13B_GSE153882_OHC", "S16_Nested_CV_summary"]],
];

for (const [filename, sheets] of files) {
  const input = await FileBlob.load(path.join(root, filename));
  const workbook = await SpreadsheetFile.importXlsx(input);
  const sheetInfo = await workbook.inspect({ kind: "sheet", include: "id,name" });
  console.log(`${filename}: ${sheetInfo.ndjson}`);

  for (const sheetName of sheets) {
    const render = await workbook.render({
      sheetName,
      range: "A1:H24",
      scale: 1,
      format: "png",
      autoCrop: "all",
    });
    const bytes = new Uint8Array(await render.arrayBuffer());
    const output = path.join(root, `${filename.replace(".xlsx", "")}_${sheetName}.png`);
    await fs.writeFile(output, bytes);
    console.log(`Rendered ${output}`);
  }
}
