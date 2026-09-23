import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const inputWorkbook = process.argv[2] ?? "LLNL_Project/results/vhh_hseq/ipi_vhh_selected_hseq.xlsx";
const payloadPath = process.argv[3] ?? "LLNL_Project/results/vhh_hseq/ipi_vhh_selected_hseq_antigen_aa_payload.json";
const outputWorkbook = process.argv[4] ?? "LLNL_Project/results/vhh_hseq/ipi_vhh_selected_hseq_antigen_aa.xlsx";

function columnLetters(index) {
  let value = index + 1;
  let letters = "";
  while (value > 0) {
    const remainder = (value - 1) % 26;
    letters = String.fromCharCode(65 + remainder) + letters;
    value = Math.floor((value - 1) / 26);
  }
  return letters;
}

const input = await FileBlob.load(inputWorkbook);
const workbook = await SpreadsheetFile.importXlsx(input);
const payload = JSON.parse(await fs.readFile(payloadPath, "utf8"));

const dataSheet = workbook.worksheets.getItem("VHH_HSEQ");
const columns = payload.columns;
const values = payload.values;
const matrix = [columns, ...values];

const startCol = 46; // Existing VHH_HSEQ table is A:AT.
const range = dataSheet.getRangeByIndexes(0, startCol, matrix.length, columns.length);
range.values = matrix;
range.format = { font: { name: "Arial", size: 10 }, wrapText: false };
dataSheet.getRangeByIndexes(0, startCol, 1, columns.length).format = {
  fill: "#8064A2",
  font: { bold: true, color: "#FFFFFF" },
};
dataSheet.getRangeByIndexes(0, startCol, matrix.length, columns.length).format.autofitColumns();
dataSheet.getRange(`${columnLetters(startCol)}:${columnLetters(startCol)}`).format.columnWidth = 80;
dataSheet.getRange(`${columnLetters(startCol + 1)}:${columnLetters(startCol + 1)}`).format.columnWidth = 80;
dataSheet.getRange(`${columnLetters(startCol + 2)}:${columnLetters(startCol + 2)}`).format.columnWidth = 60;

const lookupSheet = workbook.worksheets.add("Antigen_AA_Lookup");
lookupSheet.showGridLines = false;
const lookupMatrix = [payload.lookup_columns, ...payload.lookup_values];
const lookupLastCol = columnLetters(payload.lookup_columns.length - 1);
lookupSheet.getRangeByIndexes(0, 0, lookupMatrix.length, payload.lookup_columns.length).values = lookupMatrix;
lookupSheet.getRange(`A1:${lookupLastCol}1`).format = {
  fill: "#8064A2",
  font: { bold: true, color: "#FFFFFF" },
};
lookupSheet.getRange(`A1:${lookupLastCol}${lookupMatrix.length}`).format = {
  font: { name: "Arial", size: 10 },
  wrapText: false,
};
lookupSheet.freezePanes.freezeRows(1);
lookupSheet.getRange(`A1:${lookupLastCol}${lookupMatrix.length}`).format.autofitColumns();
lookupSheet.getRange("J:J").format.columnWidth = 80;

const summary = payload.summary;
const readme = workbook.worksheets.getItem("README");
const readmeRows = [
  ["Antigen AA Source", ""],
  ["Rows with antigen_aa", summary.rows_with_antigen_aa],
  ["Rows without antigen_aa", summary.rows_without_antigen_aa],
  ["Targets with antigen_aa", summary.targets_with_antigen_aa],
  ["Targets without antigen_aa", summary.targets_without_antigen_aa],
  ["Rows with terminal tags removed", summary.rows_with_tags_removed],
  ["Targets with terminal tags removed", summary.targets_with_tags_removed],
  ["Target AA metadata table", summary.target_aa_meta_table],
  ["AA sequence table", summary.aa_table],
  ["Tag-cleaning rule", "antigen_aa_raw preserves the construct sequence; antigen_aa trims terminal His/Strep/Avi/HA/FLAG tag tails."],
  ["Lookup rule", "Exact target match first; multiple candidates are ranked and summarized on Antigen_AA_Lookup."],
];
const readmeStart = 14;
readme.getRangeByIndexes(readmeStart - 1, 0, readmeRows.length, 2).values = readmeRows;
readme.getRange(`A${readmeStart}:B${readmeStart}`).format = {
  fill: "#8064A2",
  font: { bold: true, color: "#FFFFFF", size: 12 },
};
readme.getRange(`A${readmeStart + 1}:A${readmeStart + readmeRows.length - 1}`).format = {
  fill: "#E4DFEC",
  font: { bold: true },
};
readme.getRange(`A${readmeStart}:B${readmeStart + readmeRows.length - 1}`).format.autofitColumns();
readme.getRange("B:B").format.columnWidth = 90;

await fs.mkdir(path.dirname(outputWorkbook), { recursive: true });
const preview = await workbook.render({
  sheetName: "Antigen_AA_Lookup",
  range: "A1:L12",
  scale: 1,
  format: "png",
});
await fs.writeFile(
  outputWorkbook.replace(/\.xlsx$/i, "_lookup_preview.png"),
  new Uint8Array(await preview.arrayBuffer()),
);

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputWorkbook);
console.log(`Output workbook: ${outputWorkbook}`);
