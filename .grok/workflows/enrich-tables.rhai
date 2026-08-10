// Enrich PeerJ manuscript tables from already-computed v2 results.
// Focus: Table 1 (primary) + Table 2 (cross-tissue), optional Table 3 (per-type ranges).
// Do NOT reintroduce retired power/bootstrap/MDE inference.

let meta = #{
    name: "enrich-tables",
    description: "Enrich scReg-Eval PeerJ tables from computed v2 results without new experiments",
    when_to_use: "After figure polish, when table columns underuse already-run stats",
    phases: [
        #{ title: "Inventory", detail: "map current tables vs results JSON" },
        #{ title: "Spec", detail: "lock columns and exact numbers from JSON" },
        #{ title: "Implement", detail: "edit make_figs.R + manuscript captions" },
        #{ title: "Verify", detail: "regen tables, rebuild PDF, audit numbers" },
    ],
};

let inventory_schema = #{
    "type": "object",
    "required": ["current_tables", "unused_stats", "risks"],
    "properties": #{
        "current_tables": #{ "type": "array", "items": #{ "type": "string" } },
        "unused_stats": #{ "type": "array", "items": #{ "type": "string" } },
        "risks": #{ "type": "array", "items": #{ "type": "string" } },
        "notes": #{ "type": "string" },
    },
};

let spec_schema = #{
    "type": "object",
    "required": ["table1_plan", "table2_plan", "table3_plan", "source_files", "exact_numbers"],
    "properties": #{
        "table1_plan": #{ "type": "string" },
        "table2_plan": #{ "type": "string" },
        "table3_plan": #{ "type": "string" },
        "source_files": #{ "type": "array", "items": #{ "type": "string" } },
        "exact_numbers": #{ "type": "array", "items": #{ "type": "string" } },
        "caption_edits": #{ "type": "string" },
        "forbidden": #{ "type": "array", "items": #{ "type": "string" } },
    },
};

let implement_schema = #{
    "type": "object",
    "required": ["files_changed", "tables_written", "summary", "commands_run"],
    "properties": #{
        "files_changed": #{ "type": "array", "items": #{ "type": "string" } },
        "tables_written": #{ "type": "array", "items": #{ "type": "string" } },
        "summary": #{ "type": "string" },
        "commands_run": #{ "type": "array", "items": #{ "type": "string" } },
        "blockers": #{ "type": "array", "items": #{ "type": "string" } },
    },
};

let verify_schema = #{
    "type": "object",
    "required": ["ok", "checks", "mismatches", "pdf_pages"],
    "properties": #{
        "ok": #{ "type": "boolean" },
        "checks": #{ "type": "array", "items": #{ "type": "string" } },
        "mismatches": #{ "type": "array", "items": #{ "type": "string" } },
        "pdf_pages": #{ "type": "integer" },
        "notes": #{ "type": "string" },
    },
};

let root = if args == () { () } else { args.root };
if root == () {
    root = ".";
}
log("enrich-tables root=" + root);

// ---------------------------------------------------------------------------
phase("Inventory");
// ---------------------------------------------------------------------------
let inv_prompt = "";
inv_prompt += "You are inventorying already-computed stats vs PeerJ tables for scReg-Eval.\n";
inv_prompt += "Repo root: " + root + "\n";
inv_prompt += "Read with tools (do not guess):\n";
inv_prompt += "1) " + root + "/paper/figs/table1_primary_fixed_panel.tex\n";
inv_prompt += "2) " + root + "/paper/figs/table2_cross_tissue_observed.tex\n";
inv_prompt += "3) " + root + "/paper/make_figs.R table-writing block near end (search 'tables')\n";
inv_prompt += "4) " + root + "/paper/manuscript.tex captions for tab:primary and tab:cross\n";
inv_prompt += "5) Results: cross_tissue_additive_decomp_v2.json, cross_tissue_brain_vs_pbmc.json,\n";
inv_prompt += "   cross_tissue_atac_v2.json, cross_tissue_pbmc_vs_fibroblast.json,\n";
inv_prompt += "   fixed_panel_audit_v2.json, spec_sensitivity_v2.json,\n";
inv_prompt += "   and audit per_cell_type descriptive_summary for brain/pbmc full and non_degree.\n";
inv_prompt += "6) " + root + "/docs/LEGACY_INFERENCE_NOTE.md — list what is FORBIDDEN to put in tables.\n";
inv_prompt += "Return current_tables, unused_stats (candidates that are safe/descriptive),\n";
inv_prompt += "and risks (retired bootstrap/MDE/power, over-wide tables, PeerJ width).\n";
inv_prompt += "Prefer enriching tables over new experiments. Figures were already polished.\n";

let inv = agent(inv_prompt, #{
    label: "inventory-tables",
    capability_mode: "read-only",
    output_schema: inventory_schema,
});
if inv == () || !inv.success {
    pause("infra", "Inventory agent failed; resume after checking root path.");
}

// ---------------------------------------------------------------------------
phase("Spec");
// ---------------------------------------------------------------------------
let inv_json = json_encode(inv.output);
let spec_prompt = "";
spec_prompt += "Lock a concrete table-enrichment SPEC for the PeerJ scReg-Eval manuscript.\n";
spec_prompt += "Repo: " + root + "\n";
spec_prompt += "Inventory agent output (JSON):\n" + inv_json + "\n\n";
spec_prompt += "MANDATORY goals (user request: enrich TABLES after figure polish):\n";
spec_prompt += "A) Table 2 (cross-tissue): expand beyond Observed Spearman + Shared TFs.\n";
spec_prompt += "   Required columns (descriptive only, NO CI/bootstrap/MDE):\n";
spec_prompt += "   - Tissue pair\n";
spec_prompt += "   - Observed Spearman rho (4 decimals)\n";
spec_prompt += "   - Additive fraction explained (from decomp, as percent or 0.xx)\n";
spec_prompt += "   - Residual Spearman after additive fits\n";
spec_prompt += "   - Binary support phi\n";
spec_prompt += "   - Edge Jaccard (from pairwise cross_tissue_*.json)\n";
spec_prompt += "   - Mean per-TF-row Spearman\n";
spec_prompt += "   - Shared TFs (446)\n";
spec_prompt += "   Use display names Brain--PBMC, Brain--Fibroblast mix, PBMC--Fibroblast mix.\n";
spec_prompt += "B) Table 1 (primary full-confound): add a Support column derived from both nulls\n";
spec_prompt += "   (q_M<0.05 AND q_D<0.05 for FM rows; baseline stays outside family).\n";
spec_prompt += "   Labels: both / neither / (and if negative both) both (neg). Keep existing p/q columns.\n";
spec_prompt += "C) NEW Table 3 (optional but preferred): per-type descriptive ranges full vs non-degree\n";
spec_prompt += "   from audit per_cell_type descriptive_summary — rho min/max/median and n_rows,\n";
spec_prompt += "   brain and PBMC, full_confound and non_degree_confound. DESCRIPTIVE only.\n";
spec_prompt += "D) Update manuscript.tex captions for tab:cross, tab:primary, and add tab:pertype_ranges\n";
spec_prompt += "   if Table 3 is created. Brief in-text pointers if needed (1 sentence each).\n";
spec_prompt += "E) FORBIDDEN in any table: bootstrap CI, MDE, power, implied alpha, retired files\n";
spec_prompt += "   power_analysis_v2, cross_tissue_bootstrap_v2, stats_enhanced_v2, pertype_stats_enhanced CIs.\n\n";
spec_prompt += "Read the JSON files and compute EXACT numbers. Put each as a short string in exact_numbers.\n";
spec_prompt += "List source_files. Write table*_plan as implementable column lists.\n";
spec_prompt += "caption_edits: what to change in manuscript.tex.\n";
spec_prompt += "Do not edit files in this phase.\n";

let spec = agent(spec_prompt, #{
    label: "spec-tables",
    capability_mode: "read-only",
    output_schema: spec_schema,
});
if spec == () || !spec.success {
    pause("infra", "Spec agent failed.");
}

// Human gate: table columns change the manuscript contract.
let spec_json = json_encode(spec.output);
write_scratch_file("table_enrichment_spec.json", spec_json);
await_user("user", "Table enrichment SPEC is ready in the workflow scratch (table_enrichment_spec.json). Review columns for Table1/2/3. Resume to implement, or edit the plan first.");

// ---------------------------------------------------------------------------
phase("Implement");
// ---------------------------------------------------------------------------
let impl_prompt = "";
impl_prompt += "Implement the locked table enrichment SPEC for scReg-Eval PeerJ package.\n";
impl_prompt += "Repo root: " + root + "\n";
impl_prompt += "LOCKED SPEC JSON:\n" + spec_json + "\n\n";
impl_prompt += "Rules:\n";
impl_prompt += "1) Edit ONLY: paper/make_figs.R (table generation block), paper/manuscript.tex\n";
impl_prompt += "   (captions + optional 1-sentence in-text refs). Do NOT hand-edit generated .tex tables;\n";
impl_prompt += "   regenerate via Rscript make_figs.R.\n";
impl_prompt += "2) Table sources must be results/v2 JSON already used by make_figs (decomp, pairwise\n";
impl_prompt += "   cross_tissue_*, fixed_panel_audit primary, per_cell_type descriptive_summary).\n";
impl_prompt += "3) After R edits, run: cd " + root + "/paper && Rscript make_figs.R\n";
impl_prompt += "   Then latexmk -pdf -interaction=nonstopmode manuscript.tex\n";
impl_prompt += "4) If Table 3 is new: write figs/table3_*.tex from R, input it in manuscript.tex with\n";
impl_prompt += "   label tab:pertype_ranges near the per-type results subsection.\n";
impl_prompt += "5) Keep PeerJ width: use compact headers, 2-3 decimal places where appropriate,\n";
impl_prompt += "   fitfig wrapper like existing tables.\n";
impl_prompt += "6) Never introduce retired bootstrap/MDE numbers.\n";
impl_prompt += "7) Leave figure layout alone unless a table width forces a tiny caption tweak.\n";
impl_prompt += "Return files_changed, tables_written, commands_run, summary, blockers.\n";

let impl = agent(impl_prompt, #{
    label: "implement-tables",
    capability_mode: "all",
    output_schema: implement_schema,
});
if impl == () || !impl.success {
    pause("infra", "Implement agent failed; inspect make_figs.R / manuscript.tex.");
}
if impl.output.blockers != () && impl.output.blockers.len() > 0 {
    log("Implement reported blockers: " + json_encode(impl.output.blockers));
}

// ---------------------------------------------------------------------------
phase("Verify");
// ---------------------------------------------------------------------------
let impl_json = json_encode(impl.output);
let ver_prompt = "";
ver_prompt += "Verify table enrichment for " + root + ".\n";
ver_prompt += "Implement report:\n" + impl_json + "\n";
ver_prompt += "Spec:\n" + spec_json + "\n\n";
ver_prompt += "Do with tools:\n";
ver_prompt += "1) Read generated table tex under paper/figs/table*.tex\n";
ver_prompt += "2) Recompute expected numbers from results/v2 JSON and compare cell-by-cell\n";
ver_prompt += "3) Confirm manuscript.tex inputs the tables and captions mention descriptive-only\n";
ver_prompt += "4) Confirm NO bootstrap CI / MDE / power columns appear\n";
ver_prompt += "5) Run latexmk if PDF stale; report page count via pdfinfo or similar\n";
ver_prompt += "6) Optionally: python3 " + root + "/src/v2/build_peerj_package.py if tables changed\n";
ver_prompt += "Set ok=true only if all required Table1/Table2 enrichments match sources and PDF builds.\n";
ver_prompt += "List mismatches explicitly; empty mismatches when ok.\n";

let ver = agent(ver_prompt, #{
    label: "verify-tables",
    capability_mode: "all",
    output_schema: verify_schema,
});
if ver == () || !ver.success {
    pause("verification", "Verify agent failed.");
}
if ver.output.ok != true {
    let mm = if ver.output.mismatches == () { "unknown mismatches" } else { json_encode(ver.output.mismatches) };
    pause("verification", "Table verification failed: " + mm);
}

let report = "";
report += "# Table enrichment workflow result\n\n";
report += "## Inventory\n" + json_encode(inv.output) + "\n\n";
report += "## Spec\n" + json_encode(spec.output) + "\n\n";
report += "## Implement\n" + json_encode(impl.output) + "\n\n";
report += "## Verify\n" + json_encode(ver.output) + "\n";
let report_path = write_scratch_file("enrich-tables-report.md", report);

complete(#{
    summary: "Tables enriched and verified",
    report: report_path,
    ok: true,
    tables: impl.output.tables_written,
});
