import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const csvPath = process.argv[2] ?? "LLNL_Project/results/vhh_hseq/ipi_vhh_selected_hseq.csv";
const summaryPath = process.argv[3] ?? "LLNL_Project/results/vhh_hseq/ipi_vhh_selected_hseq_summary.json";
const outputPath = process.argv[4] ?? "LLNL_Project/results/vhh_hseq/ipi_vhh_selected_hseq.xlsx";

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

const csvText = await fs.readFile(csvPath, "utf8");
const summary = JSON.parse(await fs.readFile(summaryPath, "utf8"));
const rows = csvText.trimEnd().split(/\r?\n/);
const header = rows[0].split(",");
const rowCount = rows.length;
const lastCol = columnLetters(header.length - 1);

const workbook = await Workbook.fromCSV(csvText, { sheetName: "VHH_HSEQ" });
const dataSheet = workbook.worksheets.getItem("VHH_HSEQ");
dataSheet.showGridLines = false;
dataSheet.freezePanes.freezeRows(1);

const dataRange = dataSheet.getRange(`A1:${lastCol}${rowCount}`);
dataRange.format = {
  font: { name: "Arial", size: 10 },
  wrapText: false,
};
dataSheet.getRange(`A1:${lastCol}1`).format = {
  fill: "#1F4E79",
  font: { bold: true, color: "#FFFFFF" },
};
dataSheet.getRange(`A1:${lastCol}${rowCount}`).format.autofitColumns();
dataSheet.getRange("AT:AT").format.columnWidth = 70;

const summarySheet = workbook.worksheets.add("README");
summarySheet.showGridLines = false;
summarySheet.getRange("A1:B13").values = [
  ["VHH HSEQ Generation", ""],
  ["Source workbook", summary.input_file],
  ["Source sheet", summary.sheet],
  ["VHH scaffold library", summary.library_file],
  ["Rows processed", summary.rows],
  ["HSEQ ok rows", summary.hseq_ok_rows],
  ["HSEQ missing rows", summary.hseq_missing_rows],
  ["CDR H3 triplet rows", summary.cdr_h3_triplet_rows],
  ["CDR H3-only rows", summary.cdr_h3_only_rows],
  ["Scaffold counts", JSON.stringify(summary.scaffold_counts)],
  ["Status counts", JSON.stringify(summary.status_counts)],
  ["Triplet rule", "CDR H3 with ':' is parsed as CDR1:CDR2:CDR3."],
  ["CDR3-only rule", "CDR H3 without ':' is treated as CDR3; scaffold CDR1/CDR2 are used."],
];
summarySheet.getRange("A1:B1").format = {
  fill: "#1F4E79",
  font: { bold: true, color: "#FFFFFF", size: 12 },
};
summarySheet.getRange("A2:A13").format = {
  fill: "#D9EAF7",
  font: { bold: true },
};
summarySheet.getRange("A1:B13").format.autofitColumns();
summarySheet.getRange("B2:B13").format.columnWidth = 80;

await fs.mkdir(path.dirname(outputPath), { recursive: true });
const preview = await workbook.render({ sheetName: "README", autoCrop: "all", scale: 1, format: "png" });
await fs.writeFile(
  outputPath.replace(/\.xlsx$/i, "_readme_preview.png"),
  new Uint8Array(await preview.arrayBuffer()),
);
const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);

console.log(`Output workbook: ${outputPath}`);
