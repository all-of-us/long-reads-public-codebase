"""Tests for scripts/export_terra_configs.py. Run: python3 -m unittest discover scripts/tests"""
import json
import os
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import export_terra_configs as etc  # noqa: E402

PUBLIC = {"gcp-public-data--broad-references", "public-demo-bucket"}


def is_public(bucket):
    return bucket in PUBLIC


WS = "/workspaces/ns/ws"
LRP = "github.com/broadinstitute/long-read-pipelines/PBCCSWholeGenome"


def config(inputs, version="main", path=LRP):
    return {"namespace": "cfgns", "name": "PBCCS", "rootEntityType": "sample",
            "methodRepoMethod": {"sourceRepo": "dockstore", "methodPath": path, "methodVersion": version},
            "inputs": inputs, "outputs": {"PBCCSWholeGenome.vcf": "this.dv_vcf"}}


class FakeClient:
    def __init__(self, routes):
        self.routes = routes
        self.calls = []

    def get(self, path, allow_404=False, text=False, optional=False):
        self.calls.append(path)
        if path in self.routes:
            return self.routes[path]
        if allow_404 or optional:
            return None
        raise AssertionError(f"unexpected API call {path}")


def read(path):
    with open(path) as f:
        return f.read()


class RedactorTest(unittest.TestCase):
    WS_BUCKET = "fc-secure-0a1b2c3d-aaaa-bbbb-cccc-ddddeeeeffff"

    def setUp(self):
        self.r = etc.Redactor(is_public, {
            "ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
            "private_ref": f"gs://{self.WS_BUCKET}/refs/GRCh38_no_alt.fa",
            "tr_beds": {"itemsType": "AttributeValue", "items": ["gs://public-demo-bucket/a.bed"]},
            "sample_list": "gs://fc-secure-other/lists/1234567.txt",
        }, workspace_bucket=self.WS_BUCKET)

    def test_expressions_kept(self):
        for v in ("this.aligned_bam", "workspace.ref_fasta", "this.samples.bam"):
            self.assertEqual(self.r.value(v, "x"), v)
        self.assertEqual(self.r.events, [])

    def test_literals(self):
        self.assertEqual(self.r.value("50", "x"), "50")
        self.assertEqual(self.r.value("true", "x"), "true")
        self.assertEqual(self.r.value('"GRCh38"', "x"), '"GRCh38"')
        self.assertEqual(self.r.value("100000000", "x"), "100000000")  # numbers are never IDs here

    def test_private_bucket_masked_path_kept(self):
        out = json.loads(self.r.value(
            f'["gs://{self.WS_BUCKET}/refs/hg38.fa", "gs://fc-secure-other/train/hprc.vcf.gz", '
            f'"gs://fc-secure-other/train/hprc.bed", "gs://public-demo-bucket/x.bed"]', "x"))
        self.assertEqual(out, ["gs://<workspace-bucket>/refs/hg38.fa", "gs://<private-bucket-1>/train/hprc.vcf.gz",
                               "gs://<private-bucket-1>/train/hprc.bed", "gs://public-demo-bucket/x.bed"])
        self.assertNotIn("fc-secure", " ".join(self.r.events))

    def test_terra_bucket_without_gs_prefix_masked(self):
        cases = {
            f'"{self.WS_BUCKET}"': '"<workspace-bucket>"',
            f'"{self.WS_BUCKET}/outputs/run1"': '"<workspace-bucket>/outputs/run1"',
            f'"https://console.cloud.google.com/storage/browser/{self.WS_BUCKET}/x"':
                '"https://console.cloud.google.com/storage/browser/<workspace-bucket>/x"',
            '"fc-0a1b2c3d-1111-2222-3333-444455556666/a.txt"': '"<private-bucket-1>/a.txt"',
        }
        for raw, want in cases.items():
            self.assertEqual(self.r.value(raw, "x"), want)
        self.assertEqual(self.r.value('"gcp-fc-secure-lookalike-ok"', "x"), '"gcp-fc-secure-lookalike-ok"')

    def test_identifier_in_private_path_redacted(self):
        v = f'"gs://{self.WS_BUCKET}/samples/1234567/aligned.bam"'
        self.assertEqual(json.loads(self.r.value(v, "x")), etc.POSSIBLE_ID)

    def test_locate_reports_field_not_value(self):
        files = {"configs/a.json": {"inputs": {"W.out": '"fc-secure-1234/x"', "W.ok": "this.bam"}}}
        found = etc.locate_problems(files, is_public)
        self.assertEqual(len(found), 1)
        self.assertIn("configs/a.json: inputs.W.out", found[0])
        self.assertNotIn("1234", found[0])

    def test_ids_and_emails_redacted(self):
        self.assertEqual(json.loads(self.r.value('"sample_1234567"', "x")), etc.POSSIBLE_ID)
        self.assertEqual(json.loads(self.r.value('"1234567"', "x")), etc.POSSIBLE_ID)
        self.assertEqual(json.loads(self.r.value('"someone@broadinstitute.org"', "x")), etc.EMAIL_REDACTED)

    def test_safe_numbers_kept(self):
        for v in (f'"gs://{self.WS_BUCKET}/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"',
                  '"--sizemax 1000000 --pctseq 0.7"', '"--refdist=5000000"', '"chr1:1000000-2000000"',
                  '"chr7_75493495"'):
            self.assertNotEqual(self.r.value(v, "x"), json.dumps(etc.POSSIBLE_ID), v)
        # Still redacted: identifiers passed as arguments, or next to a safe number.
        for v in ('"--sample 1234567"', '"--sample-id=1234567"', '"--pid 1234567"', '"chr1:100-200 1234567"'):
            self.assertEqual(json.loads(self.r.value(v, "x")), etc.POSSIBLE_ID, v)
        self.assertEqual(self.r.value('"--width 1000000"', "x"), '"--width 1000000"')

    def test_text_masks_in_place(self):
        wdl = (f'String out = "gs://{self.WS_BUCKET}/out"\n'
               'meta { author: "someone@broadinstitute.org" }\nInt mem = 100000000\n')
        out = self.r.text(wdl, "x")
        self.assertIn('"gs://<workspace-bucket>/out"', out)
        self.assertIn(etc.EMAIL_REDACTED, out)
        self.assertIn("Int mem = 100000000", out)  # numbers in code are left alone

    def test_run_ids_masked(self):
        sub, wf, sub_wf = ("9f31275e-2c26-4f9b-9322-ec8a82f16081", "86f6d9d9-9a1e-401b-afc2-a7ca7dad3231",
                           "0a1b2c3d-1111-2222-3333-444455556666")
        v = f'"gs://{self.WS_BUCKET}/submissions/{sub}/Merge/{wf}/call-Concat/Sub/{sub_wf}/call-X/o.vcf.gz"'
        self.assertEqual(json.loads(self.r.value(v, "x")),
                         "gs://<workspace-bucket>/submissions/<submission-id>/Merge/<workflow-id>/"
                         "call-Concat/Sub/<workflow-id>/call-X/o.vcf.gz")
        # UUIDs outside a submission path are left alone
        self.assertEqual(etc.mask_run_ids(f"gs://public-demo-bucket/{sub}/x"), f"gs://public-demo-bucket/{sub}/x")

    def test_digits_inside_public_path_allowed(self):
        v = '"gs://gcp-public-data--broad-references/GCA_000001405.15/ref.fa"'
        self.assertEqual(self.r.value(v, "x"), v)

    def test_workspace_values(self):
        self.assertTrue(self.r.workspace_value("workspace.ref_fasta").startswith("gs://gcp-public-data"))
        self.assertEqual(self.r.workspace_value("workspace.private_ref"), "gs://<workspace-bucket>/refs/GRCh38_no_alt.fa")
        self.assertEqual(self.r.workspace_value("workspace.sample_list"), etc.POSSIBLE_ID)
        self.assertEqual(self.r.workspace_value("workspace.tr_beds"), ["gs://public-demo-bucket/a.bed"])
        self.assertIsNone(self.r.workspace_value("this.bam"))


class ExportTest(unittest.TestCase):
    def routes(self):
        current = config({"PBCCSWholeGenome.bam": "this.aligned_bam",
                          "PBCCSWholeGenome.ref": "workspace.ref_fasta",
                          "PBCCSWholeGenome.gcs_out": '"gs://fc-secure-xyz/out"',
                          "PBCCSWholeGenome.min_len": "50"})
        saved = json.loads(json.dumps(current))
        saved["inputs"]["PBCCSWholeGenome.min_len"] = "40"  # edited after the run
        old = config({"Old.bam": "this.bam"}, path="github.com/org/repo/Old")
        return {
            f"{WS}?fields=workspace.attributes,workspace.bucketName": {"workspace": {
                "bucketName": "fc-secure-xyz",
                "attributes": {"ref_fasta": "gs://gcp-public-data--broad-references/hg38/ref.fa"}}},
            f"{WS}/methodconfigs?allRepos=true": [{"namespace": "cfgns", "name": "PBCCS"}],
            f"{WS}/methodconfigs/cfgns/PBCCS": current,
            f"{WS}/submissions": [
                {"submissionId": "s1", "submissionDate": "2024-01-02T00:00:00Z", "submitter": "a@b.org",
                 "methodConfigurationNamespace": "cfgns", "methodConfigurationName": "PBCCS",
                 "workflowStatuses": {"Succeeded": 1000, "Failed": 3}},
                {"submissionId": "s0", "submissionDate": "2023-06-01T00:00:00Z",
                 "methodConfigurationNamespace": "cfgns", "methodConfigurationName": "PBCCS",
                 "workflowStatuses": {"Failed": 5}},
                {"submissionId": "s9", "submissionDate": "2023-01-01T00:00:00Z",
                 "methodConfigurationNamespace": "cfgns", "methodConfigurationName": "Deleted",
                 "workflowStatuses": {"Succeeded": 2}},
            ],
            f"{WS}/submissions/s1/configuration": saved,
            f"{WS}/submissions/s9/configuration": old,
            f"{WS}/submissions/s1": {"workflows": [
                {"workflowId": "w-failed", "status": "Failed"},
                {"workflowId": "w1", "status": "Succeeded",
                 "inputResolutions": [{"inputName": "bam", "value": "gs://fc-secure-xyz/1234567.bam"}]}]},
            f"{WS}/submissions/s1/workflows/w1?includeKey=submittedFiles": {"submittedFiles": {
                "workflowUrl": "https://raw.githubusercontent.com/org/repo/main/wdl/PBCCSWholeGenome.wdl",
                "workflow": 'version 1.0\nworkflow PBCCSWholeGenome { String out = "gs://fc-secure-xyz/o" }\n',
                "inputs": '{"PBCCSWholeGenome.bam": "gs://fc-secure-xyz/1234567.bam"}'}},
        }

    def run_export(self, routes, orch_routes=None, fetch_source=True, workspaces=("ns/ws",)):
        """Returns the exit code and the export directory of the first workspace."""
        out = tempfile.mkdtemp()
        args = SimpleNamespace(workspace=list(workspaces), out=out, api="x", orch_api="y",
                               fetch_source=fetch_source)
        rc = etc.cmd_export(args, client=FakeClient(routes), orch=FakeClient(orch_routes or {}),
                            is_public=is_public)
        return rc, os.path.join(out, etc._safe(workspaces[0]))

    def test_export(self):
        rc, out = self.run_export(self.routes())
        self.assertEqual(rc, 0)
        blob = ""
        for root, _, files in os.walk(out):
            for f in files:
                blob += open(os.path.join(root, f)).read()
        self.assertNotIn("fc-secure", blob)
        self.assertNotIn("a@b.org", blob)

        rec = json.load(open(os.path.join(out, "configs", "cfgns__PBCCS.json")))
        self.assertIn("submission s1", rec["source"])
        self.assertTrue(rec["current_differs_from_saved"])
        self.assertEqual(rec["inputs"]["PBCCSWholeGenome.min_len"], "40")  # the value actually run
        self.assertEqual(json.loads(rec["inputs"]["PBCCSWholeGenome.gcs_out"]), "gs://<workspace-bucket>/out")
        self.assertEqual(rec["workspace_values"]["workspace.ref_fasta"],
                         "gs://gcp-public-data--broad-references/hg38/ref.fa")
        self.assertEqual(rec["runs"], {"submissions": 2, "workflows_succeeded": 1000, "workflows_failed": 8,
                                       "first_run": "2023-06-01", "last_successful_run": "2024-01-02"})

        deleted = json.load(open(os.path.join(out, "configs", "cfgns__Deleted.json")))
        self.assertTrue(deleted["deleted"])
        self.assertEqual(deleted["method"]["methodPath"], "github.com/org/repo/Old")

        # The WDL that ran is kept (redacted); the resolved per-sample inputs are never written.
        src = rec["run_source"]
        self.assertTrue(src["workflow_url"].endswith("PBCCSWholeGenome.wdl"))
        wdl = read(os.path.join(out, src["workflow"]))
        self.assertIn('"gs://<workspace-bucket>/o"', wdl)
        self.assertNotIn("1234567", blob)

    def test_agora_method_source(self):
        routes = self.routes()
        agora = config({"Example.bam": "this.bam"})
        agora["methodRepoMethod"] = {"sourceRepo": "agora", "methodNamespace": "example_ns", "methodName": "ExampleMethod",
                                     "methodVersion": 2}
        routes[f"{WS}/submissions/s1/configuration"] = agora
        orch = {"/methods/example_ns/ExampleMethod/2?onlyPayload=true": "version 1.0\nworkflow ExampleMethod {}\n"}
        rc, out = self.run_export(routes, orch)
        self.assertEqual(rc, 0)
        rec = json.load(open(os.path.join(out, "configs", "cfgns__PBCCS.json")))
        self.assertIn("workflow ExampleMethod", read(os.path.join(out, rec["method_source"])))

    def test_two_workspaces(self):
        routes = self.routes()
        imp = "/workspaces/ns2/imp"
        panel = config({"Panel.vcf": '"gs://fc-secure-xyz/joint/sv.vcf.gz"',   # reads the first workspace
                        "Panel.res": '"gs://fc-secure-impbucket/resources/map.txt"'})
        routes.update({
            f"{imp}?fields=workspace.attributes,workspace.bucketName": {"workspace": {
                "bucketName": "fc-secure-impbucket", "attributes": {}}},
            f"{imp}/methodconfigs?allRepos=true": [{"namespace": "cfgns", "name": "PBCCS"}],
            f"{imp}/methodconfigs/cfgns/PBCCS": panel,
            f"{imp}/submissions": [],
        })
        rc, first = self.run_export(routes, workspaces=("ns/ws", "ns2/imp"))
        self.assertEqual(rc, 0)
        second = os.path.join(os.path.dirname(first), "ns2__imp")
        rec = json.loads(read(os.path.join(second, "configs", "cfgns__PBCCS.json")))
        self.assertEqual(json.loads(rec["inputs"]["Panel.vcf"]), "gs://<workspace-bucket:ws>/joint/sv.vcf.gz")
        self.assertEqual(json.loads(rec["inputs"]["Panel.res"]), "gs://<workspace-bucket:imp>/resources/map.txt")
        self.assertEqual(json.loads(read(os.path.join(first, "manifest.json")))["workspace"], "ns/ws")

        # organize accepts the parent directory
        records, _, workspaces = etc.load_exports([os.path.dirname(first)])
        self.assertEqual(sorted({r["_workspace"] for r in records}), ["ns/ws", "ns2/imp"])

    def test_inaccessible_workspace_fails_before_export(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_export(self.routes(), workspaces=("ns/ws", "nope/missing"))
        self.assertIn("nope/missing", str(cm.exception))

    def test_no_source(self):
        rc, out = self.run_export(self.routes(), fetch_source=False)
        self.assertEqual(rc, 0)
        self.assertFalse(os.path.exists(os.path.join(out, "sources")))

    def test_gate_blocks_leaks(self):
        with mock.patch.object(etc.Redactor, "value", lambda self, raw, where: raw):
            rc, out = self.run_export(self.routes())
        self.assertEqual(rc, 1)
        self.assertEqual(os.listdir(os.path.dirname(out)), [])

    def test_identifier_in_config_name(self):
        routes = self.routes()
        routes[f"{WS}/submissions"].append({
            "submissionId": "s5", "submissionDate": "2023-02-01", "workflowStatuses": {},
            "methodConfigurationNamespace": "cfgns", "methodConfigurationName": "rerun_1234567"})
        rc, out = self.run_export(routes)
        self.assertEqual(rc, 0)
        self.assertFalse(any("1234567" in f for f in os.listdir(os.path.join(out, "configs"))))


class CliTest(unittest.TestCase):
    def test_export_requires_workspace(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            etc.main(["export", "--out", "x"])


class OrganizeTest(unittest.TestCase):
    RAN = "version 1.0\nimport \"tasks/Utils.wdl\" as U\nworkflow PBCCSWholeGenome {\n  call U.X\n}\n"

    def setUp(self):
        self.wdl = tempfile.mkdtemp()
        for up, name, path, branch in (
                ("long-read-pipelines", "PBCCSWholeGenome", LRP, "main"),
                ("long-read-pipelines", "SexCheck", "github.com/broadinstitute/long-read-pipelines/SexCheck",
                 "sh_more_atomic_qc"),
                ("gatk-sv", "MainVcfQc", "github.com/broadinstitute/gatk-sv/MainVcfQc", "v0.26.2-beta")):
            d = os.path.join(self.wdl, up, name)
            os.makedirs(d)
            with open(os.path.join(d, "SOURCE.txt"), "w") as f:
                f.write(f"workflow: {name}\nbranch: {branch}\ndockstore: https://dockstore.org/workflows/{path}\n")
            with open(os.path.join(d, f"{name}.wdl"), "w") as f:  # repo copy: imports rewritten
                f.write(self.RAN.replace('"tasks/', '"repo/wdl/tasks/'))

    def export(self, workspace, records, sources=None):
        d = tempfile.mkdtemp()
        os.makedirs(os.path.join(d, "configs"))
        os.makedirs(os.path.join(d, "sources"))
        with open(os.path.join(d, "manifest.json"), "w") as f:
            json.dump({"workspace": workspace, "exported": "2026-09-24"}, f)
        for r in records:
            with open(os.path.join(d, "configs", etc._safe(r["config"]) + ".json"), "w") as f:
                json.dump(r, f)
        for rel, text in (sources or {}).items():
            with open(os.path.join(d, rel), "w") as f:
                f.write(text)
        return d

    @staticmethod
    def rec(name, succeeded, version, path=LRP, run_source=None):
        r = {"config": f"cfgns/{name}", "method": {"methodPath": path, "methodVersion": version},
             "runs": {"workflows_succeeded": succeeded, "last_successful_run": "2024-01-01"},
             "source": "configuration saved with submission s",
             "inputs": {"PBCCSWholeGenome.bam": "this.bam", "PBCCSWholeGenome.label": json.dumps(name),
                        "PBCCSWholeGenome.threads": "8", "PBCCSWholeGenome.unset": ""},
             "outputs": {}}
        if run_source:
            r["run_source"] = {"workflow": run_source}
        return r

    def organize(self, *exports, choose=(), exclude=None, exclude_build=()):
        args = SimpleNamespace(export_dirs=list(exports), wdl_root=self.wdl, choose=list(choose), womtool=None,
                               public_bucket=[], exclude_version=exclude, exclude_build=list(exclude_build),
                               unmatched_report=os.path.join(self.wdl, "unmatched.md"),
                               report=os.path.join(self.wdl, "report.md"))
        self.assertEqual(etc.cmd_organize(args), 0)
        return read(os.path.join(self.wdl, "report.md"))

    def chosen(self, rel):
        """Name of the configuration whose inputs were written to `rel`."""
        return json.loads(read(os.path.join(self.wdl, rel)))["PBCCSWholeGenome.label"]

    def test_choice_builds_branches_and_sources(self):
        e1 = self.export("ns/ws", [
            self.rec("PBCCS.GRCh38", 900, "main", run_source="sources/a.wdl"),
            self.rec("PBCCS.GRCh38_v2", 5000, "sl_v2_update"),   # most runs, but excluded below
            self.rec("PBCCS.GRCh38_failed", 0, "main"),
            self.rec("PBCCS.T2T", 800, "main", run_source="sources/b.wdl"),
            self.rec("other", 5, "main", path="github.com/org/repo/Unknown"),
        ], {"sources/a.wdl": self.RAN, "sources/b.wdl": self.RAN.replace("call U.X", "call U.Y")})
        e2 = self.export("ns2/imputation", [
            self.rec("qc", 7, "v0.26.2-beta", path="github.com/broadinstitute/gatk-sv/17-MainVcfQc")])
        report = self.organize(e1, e2, exclude="sl_v2")

        pb = "long-read-pipelines/PBCCSWholeGenome"
        self.assertEqual(self.chosen(f"{pb}/example_inputs.GRCh38.json"), "PBCCS.GRCh38")  # sl_v2 excluded
        self.assertEqual(self.chosen(f"{pb}/example_inputs.CHM13.json"), "PBCCS.T2T")
        self.assertIn("run WDL = repo", report)
        self.assertIn("run WDL differs from repo", report)

        self.assertEqual(self.chosen("gatk-sv/MainVcfQc/example_inputs.json"), "qc")  # 17- prefix, 2nd workspace
        self.assertIn("no configuration in the exported workspaces", report)  # SexCheck
        self.assertNotIn("github.com/org/repo/Unknown", report)  # unmatched: kept out of the repo report
        self.assertIn("github.com/org/repo/Unknown", read(os.path.join(self.wdl, "unmatched.md")))

    def test_override_and_stale_files(self):
        e = self.export("ns/ws", [self.rec("big", 1000, "main"), self.rec("small", 3, "kvg_test")])
        self.organize(e)
        self.assertEqual(self.chosen("long-read-pipelines/PBCCSWholeGenome/example_inputs.json"), "big")
        report = self.organize(e, choose=["PBCCSWholeGenome=cfgns/small"])
        self.assertEqual(self.chosen("long-read-pipelines/PBCCSWholeGenome/example_inputs.json"), "small")
        self.assertIn("branch differs from repo copy", report)

    def test_zero_success_never_chosen(self):
        e = self.export("ns/ws", [self.rec("PBCCS.GRCh38", 10, "main"), self.rec("PBCCS.T2T", 0, "main")])
        report = self.organize(e)
        pb = os.path.join(self.wdl, "long-read-pipelines/PBCCSWholeGenome")
        self.assertEqual(sorted(f for f in os.listdir(pb) if f.startswith("example_inputs")), ["example_inputs.json"])
        self.assertIn("no successful run", report)

    def test_small_build_group_not_written(self):
        e = self.export("ns/ws", [self.rec("PBCCS.GRCh38", 1000, "main"), self.rec("PBCCS.T2T.pilot", 25, "main")])
        report = self.organize(e)
        pb = os.path.join(self.wdl, "long-read-pipelines/PBCCSWholeGenome")
        self.assertEqual(sorted(f for f in os.listdir(pb) if f.startswith("example_inputs")), ["example_inputs.json"])
        self.assertIn("small run, not written", report)

    def test_exclude_build(self):
        e = self.export("ns/ws", [self.rec("PBCCS.GRCh38", 900, "main"), self.rec("PBCCS.T2T", 800, "main")])
        report = self.organize(e, exclude_build=["CHM13"])
        pb = os.path.join(self.wdl, "long-read-pipelines/PBCCSWholeGenome")
        self.assertEqual(sorted(f for f in os.listdir(pb) if f.startswith("example_inputs")), ["example_inputs.json"])
        self.assertEqual(self.chosen("long-read-pipelines/PBCCSWholeGenome/example_inputs.json"), "PBCCS.GRCh38")
        self.assertNotIn("PBCCS.T2T", report)

    def test_example_inputs_format(self):
        e = self.export("ns/ws", [self.rec("PBCCS.GRCh38", 10, "main")])
        self.organize(e)
        got = json.loads(read(os.path.join(self.wdl, "long-read-pipelines/PBCCSWholeGenome/example_inputs.json")))
        self.assertEqual(got, {"PBCCSWholeGenome.bam": "this.bam",       # Terra expression kept as a string
                               "PBCCSWholeGenome.label": "PBCCS.GRCh38",  # JSON string literal decoded
                               "PBCCSWholeGenome.threads": 8})            # JSON number decoded; unset input dropped
        self.assertFalse(any(f.endswith(".md") for f in os.listdir(os.path.join(self.wdl, "long-read-pipelines"))))

    def test_unknown_build_dropped_when_known_build_exists(self):
        test = self.rec("test_run", 3, "main")
        test["inputs"] = {"PBCCSWholeGenome.bam": '"gs://b/test.txt"'}
        e = self.export("ns/ws", [self.rec("PBCCS.GRCh38", 2, "main"), test])
        report = self.organize(e)
        pb = os.path.join(self.wdl, "long-read-pipelines/PBCCSWholeGenome")
        self.assertEqual(sorted(f for f in os.listdir(pb) if f.startswith("example_inputs")), ["example_inputs.json"])
        self.assertIn("build unknown, not written", report)

    def test_reference_build(self):
        self.assertEqual(etc.reference_build({"config": "x/PAV.T2T"}), "CHM13")
        self.assertEqual(etc.reference_build({"config": "x/y", "inputs": {"a": '"gs://b/GCA_000001405.15.fa"'}}),
                         "GRCh38")
        self.assertEqual(etc.reference_build({"config": "x/y"}), "unknown")
        self.assertEqual(etc.reference_build({"config": "x/y", "root_entity_type": "Integratedcall-t2t"}), "CHM13")
        self.assertEqual(etc.reference_build({"config": "x/y", "inputs": {
            "W.haps_vs_chm13_csv": "this.csv", "W.ref": '"gs://b/GRCh38_no_alt.fa"'}}), "GRCh38")


if __name__ == "__main__":
    unittest.main()
