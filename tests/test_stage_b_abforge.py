import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.stage_b_abforge import export_stage_b_candidates, import_stage_b_scores


VHH_SEQ = (
    "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRF"
    "TISRDNAKNTVYLQMNSLKPEDTAVYYCAARGGYYAMDYWGQGTQVTVSS"
)
ANTIGEN_SEQ = "M" + ("ACDEFGHIKLMNPQRSTVWY" * 4)


class StageBAbForgeBridgeTests(unittest.TestCase):
    def _write_demo_run(self, root: Path) -> Path:
        by_protein = root / "by_protein"
        by_protein.mkdir()
        pd.DataFrame(
            [
                {
                    "TAB_ID": "clone_001",
                    "CDR3": "AARGGYYAMDY",
                    "cdr3_aa": "AARGGYYAMDY",
                    "HSEQ": VHH_SEQ,
                    "max_freq": 0.12,
                    "mean_delphi_score": 0.81,
                    "mean_psr_score": 0.18,
                    "mean_sec_score": 0.12,
                    "library": "vhh_full",
                },
                {
                    "TAB_ID": "clone_002",
                    "CDR3": "AARTTYYAMDY",
                    "cdr3_aa": "AARTTYYAMDY",
                    "HSEQ": VHH_SEQ.replace("GGYY", "TTYY"),
                    "max_freq": 0.08,
                    "mean_delphi_score": 0.67,
                    "mean_psr_score": 0.22,
                    "mean_sec_score": 0.10,
                    "library": "vhh_full",
                },
            ]
        ).to_excel(by_protein / "TargetA_final_leads.xlsx", index=False)
        return by_protein / "TargetA_final_leads.xlsx"

    def test_export_stage_b_candidates_with_target_sequence_table(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._write_demo_run(root)
            target_table = root / "targets.csv"
            pd.DataFrame([{"target": "TargetA", "antigen_aa": ANTIGEN_SEQ}]).to_csv(target_table, index=False)

            result = export_stage_b_candidates(root, target_sequence_table=target_table)

            candidates = pd.read_csv(result.output_csv)
            ready = pd.read_csv(result.ready_csv)
            self.assertEqual(result.rows, 2)
            self.assertEqual(result.ready_rows, 2)
            self.assertEqual(len(ready), 2)
            self.assertIn("stage_b_example_id", candidates.columns)
            self.assertEqual(set(candidates["target"]), {"TargetA"})
            self.assertTrue(candidates["stage_b_ready"].all())

    def test_import_stage_b_scores_writes_enriched_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._write_demo_run(root)
            target_table = root / "targets.csv"
            pd.DataFrame([{"target": "TargetA", "antigen_aa": ANTIGEN_SEQ}]).to_csv(target_table, index=False)
            exported = export_stage_b_candidates(root, target_sequence_table=target_table)
            candidates = pd.read_csv(exported.output_csv)
            scores = candidates[["stage_b_example_id"]].copy()
            scores["interface_viability_label_pred"] = [0.91, 0.42]
            scores["binder_label_pred"] = [0.75, 0.30]
            scores["pose_success_label_pred"] = [0.88, 0.20]
            scores["iptm"] = [0.86, 0.40]
            scores_path = root / "abforge_scores.csv"
            scores.to_csv(scores_path, index=False)

            result = import_stage_b_scores(root, scores_path, target_sequence_table=target_table)

            enriched = pd.read_excel(result.combined_xlsx)
            self.assertEqual(result.rows, 2)
            self.assertEqual(result.scored_rows, 2)
            self.assertIn("abforge_structural_selection_score", enriched.columns)
            self.assertIn("abforge_selection_band", enriched.columns)
            self.assertEqual(enriched.iloc[0]["clone_id"], "clone_001")
            self.assertEqual(enriched.iloc[0]["abforge_selection_band"], "strong")
            self.assertTrue((result.output_dir / "enriched_by_protein" / "TargetA_final_leads_stage_b.xlsx").exists())


if __name__ == "__main__":
    unittest.main()
