import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.step6_ml_prediction import run_ml_prediction


class EmptyMlInputTests(unittest.TestCase):
    def test_delphi_backend_skips_empty_workbook_before_dispatch(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            by_protein = root / "by_protein"
            by_protein.mkdir()
            input_file = by_protein / "empty_target_final_leads.xlsx"
            pd.DataFrame(columns=["cdr3_aa", "vh_scaffold"]).to_excel(
                input_file, index=False
            )

            cfg = {
                "ml": {
                    "input_globs": ["by_protein/*_final_leads.xlsx"],
                    "delph": {
                        "command": ["this command must never run"],
                    },
                }
            }
            results = run_ml_prediction(root, cfg=cfg, backend="delphi")

            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["status"], "skipped_empty")
            self.assertEqual(results[0]["rows"], 0)
            self.assertFalse(
                list(by_protein.glob("*_delphi_input.xlsx")),
                "Delphi preparation must not run for an empty workbook",
            )

            status = json.loads((root / "ml_prediction_status.json").read_text())
            self.assertEqual(status[0]["status"], "skipped_empty")
            summary = pd.read_csv(root / "ml_summary.csv")
            self.assertEqual(summary.loc[0, "status"], "skipped_empty")
            self.assertEqual(summary.loc[0, "total_rows"], 0)


if __name__ == "__main__":
    unittest.main()
