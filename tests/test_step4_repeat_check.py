import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.step4_repeat_check import (
    REPEAT_RESULT_COLUMNS,
    check_repeats_for_target,
    enrich_with_repeats,
)


class EmptyRepeatCheckTests(unittest.TestCase):
    def _make_run(self, root: Path, target: str) -> None:
        (root / "by_protein").mkdir()
        (root / "cluster").mkdir()
        (root / "cluster_repeat").mkdir()

        leads = pd.DataFrame(
            columns=["cdr3_aa", "vh_scaffold", "vl_scaffold", "max_freq"]
        )
        leads.to_excel(
            root / "by_protein" / f"{target}_final_leads.xlsx", index=False
        )

        cluster = pd.DataFrame(columns=["CDR3", "heavy", "light"])
        cluster.to_excel(
            root / "cluster" / f"{target}.xlsx", sheet_name="cluster", index=False
        )

    def test_empty_target_preserves_repeat_schema_and_enriches(self):
        target = "empty_fab"
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._make_run(root, target)

            repeat_index = {"by_len": {}, "by_cdr3": {}}
            check_repeats_for_target(
                target, root, repeat_index, root / "cluster_repeat"
            )

            repeat_file = root / "cluster_repeat" / f"{target}.csv"
            repeat_columns = pd.read_csv(repeat_file).columns
            self.assertTrue(set(REPEAT_RESULT_COLUMNS).issubset(repeat_columns))

            enrich_with_repeats(root, [target])

            enriched = pd.read_excel(
                root / "by_protein" / f"{target}_final_leads.xlsx"
            )
            self.assertIn("is_repeat", enriched.columns)
            self.assertTrue(enriched.empty)
            self.assertTrue(
                (
                    root
                    / "by_protein"
                    / "final_leads_dedup_bytopfreq"
                    / f"{target}_final_leads_dedup_bytopfreq.xlsx"
                ).exists()
            )

    def test_legacy_empty_repeat_file_without_match_columns_is_safe(self):
        target = "legacy_empty_fab"
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._make_run(root, target)

            pd.DataFrame(columns=["CDR3", "heavy", "light"]).to_csv(
                root / "cluster_repeat" / f"{target}.csv", index=False
            )

            enrich_with_repeats(root, [target])

            enriched = pd.read_excel(
                root / "by_protein" / f"{target}_final_leads.xlsx"
            )
            self.assertIn("is_repeat", enriched.columns)
            self.assertTrue(enriched.empty)


if __name__ == "__main__":
    unittest.main()
