# Replication scripts: RAPID (EPJ Data Science, revision 1)

Scripts for the analyses added in the revision of *RAPID: Risk of Attribute Prediction-Induced
Disclosure* (Templ, Thees, Müller): Sections 5.1 (alternative attack scenarios, Table 3, and
risk contribution of each quasi-identifier, Figure 4),
5.2 (anonymization schemes, held-out records and Anonymeter, Table 4), the attacker-knowledge run in
Section 6.1, and the permutation-null and threshold-selection checks.

Run all scripts from the root of this repository. The UCI Adult data are downloaded on first use
into `paper/rapid/data/` (not version-controlled).

| Step | Script | Output |
|---|---|---|
| Table 3 | `B_threat_models.R` | `results/B_summary.csv` |
| Table 4, step 1 | `A1_prepare.R` | train/held-out split, CART releases |
| Table 4, step 2 | `A2_mst_anonymeter.py mst`, then `A2_mst_anonymeter.py anonymeter` | MST releases, `results/A_anonymeter.csv` |
| Table 4, step 3 | `A3_evaluate.R` | `results/A_summary.csv` |
| Figure 4 | `G_qi_contribution.R`, then `G_plot.R` | `results/G_qi_contribution_summary.csv`, `results/plot_adult_qi_contribution.pdf` |
| Null columns of Tables 3 and 4 | `E_nulls.R` | `results/E_nulls.rds` |
| Section 6.1 | `C_simulation.R` | `results/C_simulation.rds` |
| tau* (Section 7.3) | `F_threshold_select.R` | `results/F_threshold_select.rds` |
| By-class diagnostic | `D_null_check.R` | `results/D_null_check.rds` |

R: `riskutility` (>= 0.2.0; `G_qi_contribution.R` uses the paper's own loop, the same procedure is
`rapid_qi_contribution()` in the development version), `synthpop`, `ranger`, `ggplot2`. RAPID as defined in the paper uses
`cat_eval_method = "RCS_marginal"` (see `common.R`).

Python 3.11 for step 2:
`pip install smartnoise-synth anonymeter pandas scikit-learn`, then e.g.
`python -I paper/rapid/A2_mst_anonymeter.py mst`.
