"""Credential-free fixtures only; does NOT validate live LLM, UniProt or AF3."""
import argparse
from pathlib import Path

from . import library_acceptance as ev
from .prompt_workflow import Services, initialize_context, create_plan, run_plan

SMOKE_PLAN = dict(version=1, summary="OFFLINE FIXTURE: protein/PDB/AF3 preparation", clarifications=[], steps=[
    dict(id="protein", action="protein_resolve", params=dict(gene="WEE1", organism_id=9606)),
    dict(id="structures", action="pdb_search", params=dict(protein_step="protein")),
    dict(id="prediction", action="af3_prepare", params=dict(protein_step="protein", name="fixture_only",
                                                          start=2, end=8, seeds=[1], ligand_ccd=["ATP"]))])


class OfflineServices(Services):
    mode = "offline_synthetic_fixtures_not_live_evidence"

    @staticmethod
    def fetch_json(url):
        if "/search?" in url: return dict(results=[dict(primaryAccession="P30291")])
        if "/chemcomp/ATP" in url: return dict(chem_comp=dict(id="ATP"), fixture=True)
        if url.endswith("/P30291.json"):
            return dict(primaryAccession="P30291", entryType="UniProtKB reviewed (Swiss-Prot)",
                        organism=dict(taxonId=9606, scientificName="fixture species"),
                        genes=[dict(geneName=dict(value="WEE1"))],
                        sequence=dict(value="ACDEFGHIKLMNPQRSTVWY", length=20), fixture=True)
        raise AssertionError("Unexpected offline request")

    @staticmethod
    def pdb_search(accession, **kwargs):
        return {"fixture": True, "accession": accession}, [{"pdb_id": "8BJU", "fixture": True}]

    @staticmethod
    def run_command(*_):
        raise AssertionError("Offline preparation smoke must not execute compute")


def run_smoke(output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    context = initialize_context(output / "workspace", "smoke", "E038 offline smoke")
    fixture = output / "fixture-plan.json"; ev.write(fixture, SMOKE_PLAN)
    plan = create_plan(Path(context["db"]), context["user_id"], context["project_id"],
                       "Offline synthetic fixture only", response_file=fixture)
    services = OfflineServices()
    report, path = run_plan(Path(context["db"]), context["user_id"], context["project_id"], plan, services=services)
    if report["status"] != "complete": raise ValueError(f"Offline workflow failed; inspect {path}")
    prepared = ev.read(path.parent / "prediction/af3-input.json")
    if prepared["sequences"][0]["protein"]["sequence"] != "CDEFGHI": raise ValueError("Construct mismatch")
    repeat, _ = run_plan(Path(context["db"]), context["user_id"], context["project_id"], plan, services=services)
    if repeat != report: raise ValueError("Resume result changed")
    summary = dict(status="passed", mode=services.mode, project=context["project_id"], report=str(path),
                   live_llm="not_run", live_protein_api="not_run", af3_inference="not_run", real_3d_search="not_run",
                   verified=["Project scope", "structured plan", "typed dependencies", "sequence construct",
                             "CCD input", "provenance", "hash-checked resume"])
    ev.write(output / "report.json", summary)
    print(f"Offline smoke passed: {output / 'report.json'}")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run_smoke(parser.parse_args().output)


if __name__ == "__main__": main()
