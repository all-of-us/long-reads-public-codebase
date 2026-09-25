#!/usr/bin/env python3
"""Export the Terra workflow configurations used for Phase 1, for publication.

Two subcommands, run in two places:

  export    Run on a Terra VM. For each Phase 1 workspace, reads every workflow
            configuration and the submission history through the Terra API,
            plus the WDL that each configuration's latest successful run
            actually submitted, and the WDL of any Terra Methods Repository
            method. Redacts anything not public and writes a directory for
            review. It never reads per-sample data tables, submitter
            identities, or resolved per-sample inputs.

  organize  Run locally from the repository root on one or more reviewed
            exports. Matches each configuration to a workflow under
            phase1/terra/wdl/, picks one configuration per workflow and
            reference build, and writes its inputs next to the workflow as
            example_inputs.json. A review table of the choices and the
            alternatives is written next to the export, not into the repository.

Typical use:

  # on a Terra Jupyter VM with access to the workspaces to export
  # (Python 3.8+, standard library only; CI tests 3.10)
  python3 export_terra_configs.py export --workspace NAMESPACE/NAME [--workspace ...] \
      --out terra_config_export
  #   ...review terra_config_export/<workspace>/, then download terra_config_export/...

  # locally, from the repository root
  python3 scripts/export_terra_configs.py organize terra_config_export --exclude-version sl_v2

Redaction. Terra stores each input as a string: an expression such as
`this.aligned_bam` or `workspace.ref_fasta`, or a JSON literal. Expressions are
kept, because they carry no participant data and show how outputs feed inputs.
In literal values, resolved `workspace.*` values, and fetched WDL source, the
name of every bucket that is not anonymously readable is replaced by a label,
but the rest of the path is kept: gs://fc-secure-<uuid>/refs/GRCh38.fa becomes
gs://<workspace-bucket>/refs/GRCh38.fa. Other private buckets become
<private-bucket-1>, <private-bucket-2>, ..., numbered consistently across the
export. Terra submission and workflow IDs in paths (gs://.../submissions/<uuid>/Workflow/
<uuid>/call-...) become <submission-id> and <workflow-id>. Email addresses are removed. An input value is replaced entirely if it
holds a run of 7+ digits (a possible participant ID), except in public paths,
genome accessions (GCA_000001405.15), genomic coordinates (chr1:1000000-2000000),
and numeric command-line arguments (--sizemax 1000000, but not arguments named
like sample, name, or id). WDL source is checked for buckets and emails only.
As a final gate, nothing is written if the output still contains a private
bucket name or an email address.
"""
import argparse
import datetime
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from typing import Callable, Dict, Iterable, List, Optional

DEFAULT_API = "https://rawls.dsde-prod.broadinstitute.org/api"
DEFAULT_ORCH_API = "https://api.firecloud.org/api"
KNOWN_PUBLIC_BUCKETS = {
    "gcp-public-data--broad-references",
    "genomics-public-data",
    "broad-public-datasets",
    "gcp-public-data--gnomad",
}

GS_URI = re.compile(r"gs://([a-z0-9][a-z0-9._-]*[a-z0-9])(/[^\s\"',\]]*)?")
HTTPS_GCS = re.compile(r"https://storage\.(?:googleapis|cloud\.google)\.com/([a-z0-9][a-z0-9._-]*[a-z0-9])/")
# Terra workspace buckets are fc-<uuid> or fc-secure-<uuid>. Catch them without a gs:// prefix too,
# for example a bare bucket name or a console.cloud.google.com/storage/browser/ URL.
TERRA_BUCKET = re.compile(r"(?<![A-Za-z0-9_-])(fc-secure-[a-z0-9-]+|fc-[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12})"
                          r"(/[^\s\"',\]]*)?")
UUID = r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
LONG_DIGITS = re.compile(r"(?<![0-9A-Za-z])\d{7,}(?![0-9A-Za-z])")
# Numbers that are not identifiers: genome accessions, genomic coordinates, numeric CLI arguments.
SAFE_NUMERIC = re.compile(r"GC[AF]_\d+(?:\.\d+)?"
                          r"|chr[0-9XYMT]{1,2}[:_]\d+(?:-\d+)?"
                          r"|(?<![\w.])--?(?![\w-]*(?:sample|name|person|participant|subject|id\b))"
                          r"[A-Za-z][\w-]*(?:=|\s+)\d+(?:\.\d+)?")

WORKSPACE_BUCKET = "<workspace-bucket>"
POSSIBLE_ID = "<redacted: possible participant identifier>"
EMAIL_REDACTED = "<redacted: email>"


# --------------------------------------------------------------------------- #
# Redaction
# --------------------------------------------------------------------------- #

def mask_run_ids(s: str) -> str:
    """Replace Terra submission and workflow IDs in execution paths with placeholders."""
    if "/submissions/" not in s:
        return s
    s = re.sub(r"(/submissions/)" + UUID, r"\1<submission-id>", s)
    return re.sub(r"(?<=/)" + UUID + r"(?=/)", "<workflow-id>", s)


def _map_strings(v, f: Callable[[str], str]):
    if isinstance(v, str):
        return f(v)
    if isinstance(v, list):
        return [_map_strings(x, f) for x in v]
    if isinstance(v, dict):
        return {k: _map_strings(x, f) for k, x in v.items()}
    return v


class Redactor:
    """Redacts Terra input/output strings and WDL text. `is_public(bucket)` decides bucket visibility."""

    def __init__(self, is_public: Callable[[str], bool], workspace_attrs: Dict[str, object],
                 workspace_bucket: Optional[str] = None, workspace_label: str = WORKSPACE_BUCKET,
                 labels: Optional[Dict[str, str]] = None):
        """`labels` may be shared between Redactors so a bucket gets the same label in every workspace."""
        self.is_public = is_public
        self.workspace_attrs = workspace_attrs
        self.labels: Dict[str, str] = labels if labels is not None else {}
        if workspace_bucket:
            self.labels[workspace_bucket] = workspace_label
        self.events: List[str] = []

    def label(self, bucket: str) -> str:
        if bucket not in self.labels:
            n = sum(1 for v in self.labels.values() if v.startswith("<private-bucket-")) + 1
            self.labels[bucket] = f"<private-bucket-{n}>"
        return self.labels[bucket]

    def _mask_buckets(self, s: str, where: str) -> str:
        masked = []

        def gs(m: "re.Match") -> str:
            if self.is_public(m.group(1)):
                return m.group(0)
            masked.append(1)
            return f"gs://{self.label(m.group(1))}{m.group(2) or ''}"

        def https(m: "re.Match") -> str:
            if self.is_public(m.group(1)):
                return m.group(0)
            masked.append(1)
            return f"https://storage.googleapis.com/{self.label(m.group(1))}/"

        def bare(m: "re.Match") -> str:
            masked.append(1)
            return f"{self.label(m.group(1))}{m.group(2) or ''}"

        s = TERRA_BUCKET.sub(bare, HTTPS_GCS.sub(https, GS_URI.sub(gs, s)))
        if masked:
            self.events.append(f"{where}: bucket name masked")
        return s

    def _string(self, s: str, where: str) -> str:
        if EMAIL.search(s):
            self.events.append(f"{where}: email address")
            return EMAIL_REDACTED
        s = mask_run_ids(self._mask_buckets(s, where))
        # Only public gs:// paths are left for GS_URI to match; digit runs in them are fine. Digits
        # anywhere else, including in masked private paths, are checked after removing safe numbers.
        if LONG_DIGITS.search(SAFE_NUMERIC.sub("", GS_URI.sub("", s))):
            self.events.append(f"{where}: 7+ digit identifier")
            return POSSIBLE_ID
        return s

    def _json(self, v, where: str):
        if isinstance(v, str):
            return self._string(v, where)
        if isinstance(v, list):
            return [self._json(x, where) for x in v]
        if isinstance(v, dict):
            return {k: self._json(x, where) for k, x in v.items()}
        return v  # numbers, booleans, null

    def value(self, raw: str, where: str) -> str:
        """Redact one Terra input or output string; returns a string in the same form."""
        stripped = raw.strip()
        if stripped == "" or re.fullmatch(r"(this|workspace)(\.[A-Za-z0-9_\-]+)+", stripped):
            return raw
        try:
            parsed = json.loads(stripped)
        except ValueError:
            return self._string(raw, where)  # not JSON: treat as a bare string
        return json.dumps(self._json(parsed, where))

    def text(self, s: str, where: str) -> str:
        """Redact free text such as WDL source: mask buckets and remove emails in place."""
        if EMAIL.search(s):
            self.events.append(f"{where}: email address")
            s = EMAIL.sub(EMAIL_REDACTED, s)
        return self._mask_buckets(s, where)

    def workspace_value(self, expr: str) -> Optional[object]:
        """Resolved, redacted value of `workspace.X` when it is a path or list of paths."""
        m = re.fullmatch(r"workspace\.([A-Za-z0-9_\-]+)", expr.strip())
        if not m or m.group(1) not in self.workspace_attrs:
            return None
        val = self.workspace_attrs[m.group(1)]
        if isinstance(val, dict) and "items" in val:  # Terra list attribute
            val = val["items"]
        if not isinstance(val, (str, list)) or "gs://" not in json.dumps(val):
            return None
        return self._json(val, expr.strip())


def final_gate(text: str, is_public: Callable[[str], bool]) -> List[str]:
    """Problems that must block writing."""
    problems = []
    for m in list(GS_URI.finditer(text)) + list(HTTPS_GCS.finditer(text)):
        if not is_public(m.group(1)):
            problems.append("private bucket path")
    if "fc-secure" in text or TERRA_BUCKET.search(text):
        problems.append("Terra workspace bucket name")
    if EMAIL.search(text):
        problems.append("email address")
    return sorted(set(problems))


def locate_problems(files: Dict[str, object], is_public: Callable[[str], bool]) -> List[str]:
    """Where each gate problem is, as 'file: field (problem)'. Never includes the value itself."""
    found = []

    def walk(obj, where):
        if isinstance(obj, dict):
            for k, v in obj.items():
                key_problems = final_gate(k, is_public)
                if key_problems:
                    found.append(f"{where}: a key name ({', '.join(key_problems)})")
                walk(v, f"{where}.{k}" if where.count(":") else f"{where}: {k}")
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                walk(v, f"{where}[{i}]")
        elif isinstance(obj, str):
            problems = final_gate(obj, is_public)
            if problems:
                found.append(f"{where} ({', '.join(problems)})")

    for rel, content in files.items():
        walk(content, rel)
    return found


def report_problems(header: str, files: Dict[str, object], is_public: Callable[[str], bool]) -> bool:
    """Print located problems to stderr; True if any were found."""
    blob = "\n".join(v if isinstance(v, str) else json.dumps(v) for v in files.values())
    if not final_gate(blob, is_public):
        return False
    print(header, file=sys.stderr)
    for p in locate_problems(files, is_public) or ["(could not locate; please report)"]:
        print(f"  - {p}", file=sys.stderr)
    return True


# --------------------------------------------------------------------------- #
# Terra API
# --------------------------------------------------------------------------- #

class TerraClient:
    def __init__(self, api: str, token: str):
        self.api = api.rstrip("/")
        self.token = token

    def get(self, path: str, allow_404: bool = False, text: bool = False, optional: bool = False):
        """GET a path. allow_404: None on 404/405. optional: also None on 403 (no access)."""
        req = urllib.request.Request(self.api + path, headers={
            "Authorization": f"Bearer {self.token}", "Accept": "text/plain" if text else "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                body = r.read().decode("utf-8")
                return body if text else json.loads(body)
        except urllib.error.HTTPError as e:
            if (allow_404 or optional) and e.code in (404, 405) or optional and e.code == 403:
                return None
            raise SystemExit(f"Terra API {e.code} for {path}: {e.read()[:300]!r}")


def workspace_path(ns: str, name: str) -> str:
    return f"/workspaces/{urllib.parse.quote(ns)}/{urllib.parse.quote(name)}"


def anonymous_bucket_check() -> Callable[[str], bool]:
    """True if the bucket's object listing is readable without credentials."""
    cache: Dict[str, bool] = {}

    def is_public(bucket: str) -> bool:
        if bucket in KNOWN_PUBLIC_BUCKETS:
            return True
        if bucket.startswith("fc-"):
            return False  # Terra workspace buckets
        if bucket not in cache:
            url = f"https://storage.googleapis.com/storage/v1/b/{bucket}/o?maxResults=1"
            try:
                with urllib.request.urlopen(url, timeout=30) as r:
                    cache[bucket] = r.status == 200
            except (urllib.error.URLError, OSError):
                cache[bucket] = False
        return cache[bucket]

    return is_public


def access_token() -> str:
    if os.environ.get("TERRA_TOKEN"):
        return os.environ["TERRA_TOKEN"]
    try:
        return subprocess.check_output(["gcloud", "auth", "print-access-token"], text=True).strip()
    except (OSError, subprocess.CalledProcessError) as e:
        raise SystemExit(f"Could not get an access token from gcloud ({e}). Set TERRA_TOKEN instead.")


# --------------------------------------------------------------------------- #
# export
# --------------------------------------------------------------------------- #

class SourceStore:
    """Deduplicated WDL texts, written as sources/<hash>.wdl."""

    def __init__(self):
        self.texts: Dict[str, str] = {}

    def add(self, text: str) -> str:
        rel = f"sources/{hashlib.sha256(text.encode()).hexdigest()[:16]}.wdl"
        self.texts[rel] = text
        return rel


def _date(s: Optional[str]) -> str:
    return (s or "")[:10]


def _succeeded(s: dict) -> int:
    return (s.get("workflowStatuses") or {}).get("Succeeded", 0)


def summarize_runs(subs: List[dict]) -> dict:
    ok = [s for s in subs if _succeeded(s) > 0]
    return {
        "submissions": len(subs),
        "workflows_succeeded": sum(_succeeded(s) for s in subs),
        "workflows_failed": sum((s.get("workflowStatuses") or {}).get("Failed", 0) for s in subs),
        "first_run": min((_date(s.get("submissionDate")) for s in subs), default=""),
        "last_successful_run": max((_date(s.get("submissionDate")) for s in ok), default=""),
    }


def latest_successful(subs: List[dict]) -> Optional[dict]:
    ok = sorted((s for s in subs if _succeeded(s) > 0), key=lambda s: s.get("submissionDate", ""), reverse=True)
    return ok[0] if ok else None


def export_config(client, ws: str, key: str, current: Optional[dict], subs: List[dict], red: Redactor) -> dict:
    """Prefer the configuration saved with the latest successful submission."""
    latest = latest_successful(subs)
    snapshot, source = None, None
    if latest:
        snapshot = client.get(f"{ws}/submissions/{latest['submissionId']}/configuration", allow_404=True)
        if snapshot:
            source = f"configuration saved with submission {latest['submissionId']} ({_date(latest.get('submissionDate'))})"
    cfg = snapshot or current
    if cfg is None:
        return {"config": key, "error": "configuration deleted and no submission snapshot available",
                "runs": summarize_runs(subs)}
    if source is None:
        source = "current configuration (no saved submission configuration was available)"

    differs = None
    if snapshot and current:
        differs = (snapshot.get("inputs") != current.get("inputs")
                   or snapshot.get("outputs") != current.get("outputs"))

    method = cfg.get("methodRepoMethod") or {}
    inputs = {k: red.value(v, f"{key} input {k}") for k, v in sorted((cfg.get("inputs") or {}).items())}
    outputs = {k: red.value(v, f"{key} output {k}") for k, v in sorted((cfg.get("outputs") or {}).items())}
    ws_values = {}
    for v in list((cfg.get("inputs") or {}).values()):
        resolved = red.workspace_value(v)
        if resolved is not None:
            ws_values[v.strip()] = resolved
    return {
        "config": key,
        "deleted": current is None,
        "method": {k: method.get(k) for k in ("sourceRepo", "methodPath", "methodVersion",
                                               "methodNamespace", "methodName") if method.get(k)},
        "root_entity_type": cfg.get("rootEntityType"),
        "source": source,
        "current_differs_from_saved": differs,
        "runs": summarize_runs(subs),
        "inputs": inputs,
        "outputs": outputs,
        "workspace_values": dict(sorted(ws_values.items())),
    }


def fetch_run_source(client, ws: str, sid: str, red: Redactor, store: SourceStore, where: str) -> Optional[dict]:
    """WDL submitted by one successful workflow of a submission. Reads only `submittedFiles`, and from it
    only the workflow source, its URL, and imports; never the resolved inputs."""
    detail = client.get(f"{ws}/submissions/{sid}", optional=True) or {}
    wid = next((w.get("workflowId") for w in detail.get("workflows", [])
                if w.get("status") == "Succeeded" and w.get("workflowId")), None)
    del detail  # holds per-sample inputResolutions; drop it immediately
    if not wid:
        return None
    meta = client.get(f"{ws}/submissions/{sid}/workflows/{wid}?includeKey=submittedFiles", optional=True) or {}
    sf = meta.get("submittedFiles") or {}
    out = {}
    if isinstance(sf.get("workflowUrl"), str) and sf["workflowUrl"]:
        out["workflow_url"] = red.text(sf["workflowUrl"], f"{where} workflow URL")
    if isinstance(sf.get("workflow"), str) and sf["workflow"].strip():
        out["workflow"] = store.add(red.text(sf["workflow"], f"{where} workflow source"))
    if isinstance(sf.get("imports"), dict):
        out["imports"] = {red.text(k, f"{where} import"): store.add(red.text(v, f"{where} import"))
                          for k, v in sf["imports"].items() if isinstance(v, str)}
    return out or None


def fetch_agora_source(orch, method: dict, red: Redactor, store: SourceStore, cache: Dict) -> Optional[str]:
    ns, name, ver = method.get("methodNamespace"), method.get("methodName"), method.get("methodVersion")
    if method.get("sourceRepo") != "agora" or not (ns and name and ver):
        return None
    k = (ns, name, str(ver))
    if k not in cache:
        q = urllib.parse.quote
        body = orch.get(f"/methods/{q(ns)}/{q(name)}/{q(str(ver))}?onlyPayload=true", text=True, optional=True)
        if body is not None:
            try:
                parsed = json.loads(body)
                body = parsed.get("payload") if isinstance(parsed, dict) else parsed
            except ValueError:
                pass
        cache[k] = store.add(red.text(body, f"method {ns}/{name}/{ver}")) if isinstance(body, str) else None
    return cache[k]


def export_workspace(client, orch, ns: str, name: str, red: Redactor, fetch_source: bool) -> Dict[str, object]:
    """All files for one workspace's export, keyed by path relative to its export directory."""
    ws = workspace_path(ns, name)
    store = SourceStore()
    agora_cache: Dict = {}
    configs = {f"{c['namespace']}/{c['name']}": c
               for c in client.get(f"{ws}/methodconfigs?allRepos=true") or []}
    subs_by_cfg: Dict[str, List[dict]] = {}
    for s in client.get(f"{ws}/submissions") or []:
        key = f"{s.get('methodConfigurationNamespace')}/{s.get('methodConfigurationName')}"
        subs_by_cfg.setdefault(key, []).append(s)

    keys = sorted(set(configs) | set(subs_by_cfg))
    records = []
    for i, key in enumerate(keys):
        if i % 100 == 0:
            print(f"  {name}: {i}/{len(keys)} configurations...", file=sys.stderr)
        current = None
        if key in configs:
            cns, cname = key.split("/", 1)
            current = client.get(f"{ws}/methodconfigs/{urllib.parse.quote(cns)}/{urllib.parse.quote(cname)}",
                                 allow_404=True)
        subs = subs_by_cfg.get(key, [])
        record = export_config(client, ws, key, current, subs, red)
        if fetch_source and "error" not in record:
            latest = latest_successful(subs)
            if latest:
                run_source = fetch_run_source(client, ws, latest["submissionId"], red, store, key)
                if run_source:
                    record["run_source"] = run_source
            method_source = fetch_agora_source(orch, record["method"], red, store, agora_cache)
            if method_source:
                record["method_source"] = method_source
        if EMAIL.search(key) or LONG_DIGITS.search(key):
            record["config"] = f"redacted-config-{i}"
            red.events.append(f"configuration name {i}: email or 7+ digit identifier")
        records.append(record)

    manifest = {
        "workspace": f"{ns}/{name}",
        "exported": datetime.date.today().isoformat(),
        "configurations": len(records),
        "configurations_with_run_source": sum(1 for r in records if "run_source" in r),
        "wdl_sources": len(store.texts),
        "redactions": red.events,
    }
    files: Dict[str, object] = {f"configs/{_safe(r['config'])}.json": r for r in records}
    files.update(store.texts)
    files["manifest.json"] = manifest
    return files


def cmd_export(args, client=None, orch=None, is_public=None) -> int:
    workspaces = []
    for w in args.workspace:
        if w.count("/") != 1:
            raise SystemExit(f"--workspace must be NAMESPACE/NAME, got {w!r}")
        workspaces.append(tuple(w.split("/")))
    if client is None or (orch is None and args.fetch_source):
        token = access_token()
        client = client or TerraClient(args.api, token)
        orch = orch or TerraClient(args.orch_api, token)
    is_public = is_public or anonymous_bucket_check()

    # Check access to every workspace before the long part.
    infos = {}
    for ns, name in workspaces:
        info = client.get(f"{workspace_path(ns, name)}?fields=workspace.attributes,workspace.bucketName",
                          optional=True)
        if info is None:
            raise SystemExit(f"No access to workspace {ns}/{name} (or it does not exist). Nothing was exported.")
        infos[(ns, name)] = info.get("workspace", {})

    labels: Dict[str, str] = {}  # shared, so a bucket keeps its label across workspaces
    for (ns, name), info in infos.items():
        if info.get("bucketName"):
            labels[info["bucketName"]] = f"<workspace-bucket:{name}>" if len(workspaces) > 1 else WORKSPACE_BUCKET

    exports: Dict[str, Dict[str, object]] = {}
    for ns, name in workspaces:
        red = Redactor(is_public, infos[(ns, name)].get("attributes", {}), labels=labels)
        files = export_workspace(client, orch, ns, name, red, args.fetch_source)
        exports[_safe(f"{ns}/{name}")] = files
    for files in exports.values():
        files["manifest.json"]["bucket_labels"] = sorted(set(labels.values()))

    everything = {f"{sub}/{rel}": content for sub, files in exports.items() for rel, content in files.items()}
    if report_problems("Refusing to write: the export still contains data that must not leave the workspace "
                       "(values are not shown):", everything, is_public):
        print("Nothing was written. Send these locations to the maintainer so the redaction rules can be "
              "extended.", file=sys.stderr)
        return 1

    for rel, content in everything.items():
        path = os.path.join(args.out, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            if isinstance(content, str):
                f.write(content)
            else:
                json.dump(content, f, indent=2)
                f.write("\n")
    for sub, files in exports.items():
        m = files["manifest.json"]
        print(f"{m['workspace']}: {m['configurations']} configurations, {m['wdl_sources']} WDL sources, "
              f"{len(m['redactions'])} values redacted -> {os.path.join(args.out, sub)}/")
    print("Review the files (especially each manifest.json 'redactions' and sources/) before copying them out.")
    return 0


def _safe(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "__", key)


# --------------------------------------------------------------------------- #
# organize
# --------------------------------------------------------------------------- #

def _read(path: str) -> str:
    with open(path) as f:
        return f.read()


def load_repo_workflows(wdl_root: str) -> List[dict]:
    out = []
    for src in sorted(_glob(wdl_root, "SOURCE.txt")):
        d = os.path.dirname(src)
        kv = dict(l.split(": ", 1) for l in _read(src).splitlines() if ": " in l)
        dockstore = kv.get("dockstore", "")
        out.append({
            "dir": d,
            "name": os.path.basename(d),
            "names": {n.lower() for n in (os.path.basename(d), kv.get("workflow"), kv.get("lookup_name")) if n},
            "dockstore_path": dockstore.split("/workflows/", 1)[1] if "/workflows/" in dockstore else None,
            "branch": kv.get("branch"),
        })
    return out


def _glob(root: str, filename: str) -> Iterable[str]:
    for dirpath, _, files in os.walk(root):
        if filename in files and "/repo" not in dirpath:
            yield os.path.join(dirpath, filename)


def _last_component(path: str) -> str:
    # Dockstore paths can carry an ordering prefix, e.g. github.com/broadinstitute/gatk-sv/17-MainVcfQc
    return re.sub(r"^\d+-", "", path.rsplit("/", 1)[-1]).lower()


def match(record: dict, workflows: List[dict]) -> Optional[dict]:
    method = record.get("method", {})
    path = (method.get("methodPath") or "").strip("/")
    for w in workflows:
        if path and w["dockstore_path"] and path.lower() == w["dockstore_path"].lower():
            return w
    if path:
        repo_prefix = path.rsplit("/", 1)[0].lower()
        for w in workflows:
            dp = (w["dockstore_path"] or "").lower()
            if dp and dp.rsplit("/", 1)[0] == repo_prefix and _last_component(dp) == _last_component(path):
                return w
    last = _last_component(path) if path else (method.get("methodName") or "").lower()
    hits = [w for w in workflows if last and last in w["names"]]
    return hits[0] if len(hits) == 1 else None


def reference_build(r: dict) -> str:
    """GRCh38, CHM13, or unknown: from the configuration name, then the table it ran on, then input values
    (values only; input names such as `haps_vs_chm13_csv` say nothing about the build)."""
    def classify(text: str) -> Optional[str]:
        t2t = re.search(r"t2t|chm13", text, re.I)
        g38 = re.search(r"grch38|hg38|gca_000001405", text, re.I)
        if t2t and not g38:
            return "CHM13"
        if g38 and not t2t:
            return "GRCh38"
        return None
    return (classify(r.get("config", ""))
            or classify(r.get("root_entity_type") or "")
            or classify(json.dumps([list(r.get("inputs", {}).values()), r.get("workspace_values", {})]))
            or "unknown")


def rank_key(r: dict):
    runs = r.get("runs", {})
    return (runs.get("workflows_succeeded", 0), runs.get("last_successful_run", ""))


def _wdl_body(text: str) -> List[str]:
    """WDL lines for comparison: imports (rewritten in the repository copy) and blank lines ignored."""
    return [l.strip() for l in text.splitlines() if l.strip() and not l.strip().startswith("import ")]


def compare_run_source(r: dict, entry: str) -> str:
    rel = (r.get("run_source") or {}).get("workflow")
    if not rel:
        url = (r.get("run_source") or {}).get("workflow_url")
        return "run WDL not recorded" + (f" (ran from {url})" if url else "")
    ran, repo = _wdl_body(_read(os.path.join(r["_export"], rel))), _wdl_body(_read(entry))
    if ran == repo:
        return "run WDL = repo"
    return f"run WDL differs from repo ({len(set(ran) ^ set(repo))} differing lines)"


def womtool_inputs(entry: str, womtool: Optional[str]) -> Optional[set]:
    if not womtool:
        return None
    cmd = ["java", "-jar", womtool] if womtool.endswith(".jar") else [womtool]
    try:
        out = subprocess.check_output(cmd + ["inputs", entry], text=True, stderr=subprocess.DEVNULL)
        return set(json.loads(out))
    except (OSError, subprocess.CalledProcessError, ValueError):
        return None


def _export_dirs(paths: List[str]) -> List[str]:
    """Each path is an export directory, or a parent holding one export directory per workspace."""
    out = []
    for p in paths:
        if os.path.exists(os.path.join(p, "manifest.json")):
            out.append(p)
        else:
            subs = sorted(os.path.dirname(m) for m in glob.glob(os.path.join(p, "*", "manifest.json")))
            if not subs:
                raise SystemExit(f"No manifest.json in {p} or its subdirectories.")
            out.extend(subs)
    return out


def load_exports(export_dirs: List[str]):
    records, files, workspaces = [], {}, []
    for d in _export_dirs(export_dirs):
        manifest = json.loads(_read(os.path.join(d, "manifest.json")))
        workspaces.append(f"`{manifest['workspace']}` (exported {manifest['exported']})")
        for path in sorted(glob.glob(os.path.join(d, "configs", "*.json"))):
            r = json.loads(_read(path))
            files[os.path.relpath(path, d)] = r
            r["_workspace"], r["_export"] = manifest["workspace"], d
            records.append(r)
        for path in sorted(glob.glob(os.path.join(d, "sources", "*.wdl"))):
            files[os.path.relpath(path, d)] = _read(path)
    return records, files, workspaces


def example_inputs(r: dict) -> dict:
    """The configuration's inputs in `Workflow.input` form: JSON literals become JSON values, Terra
    expressions (this.<column>, workspace.<name>) stay as strings, and unset inputs are dropped."""
    out = {}
    for k, v in sorted(r.get("inputs", {}).items()):
        if not v.strip():
            continue
        try:
            out[k] = json.loads(v)
        except ValueError:
            out[k] = v
    return _map_strings(out, mask_run_ids)


def cmd_organize(args) -> int:
    records, files, workspaces = load_exports(args.export_dirs)

    # Re-check the reviewed exports before anything enters the repository.
    check = anonymous_bucket_check()
    is_public = lambda b: b in args.public_bucket or check(b)  # noqa: E731
    if report_problems("Export contains data that must not enter the repository:", files, is_public):
        print("Re-run export, or pass --public-bucket NAME for a bucket you have confirmed is public.",
              file=sys.stderr)
        return 1

    workflows = load_repo_workflows(args.wdl_root)
    overrides = set(o.split("=", 1)[1] for o in args.choose)
    by_wf: Dict[str, Dict[str, List[dict]]] = {}
    unmatched = []
    excluded = 0
    for r in records:
        version = str(r.get("method", {}).get("methodVersion", ""))
        if args.exclude_version and re.search(args.exclude_version, version):
            excluded += 1
            continue
        w = match(r, workflows) if "error" not in r else None
        if w and reference_build(r) in args.exclude_build:
            excluded += 1
            continue
        if w:
            by_wf.setdefault(w["name"], {}).setdefault(reference_build(r), []).append(r)
        else:
            unmatched.append(r)

    lines = [
        "# Phase 1 Terra run configurations",
        "",
        "Generated by `scripts/export_terra_configs.py organize` from exports of "
        + ", ".join(workspaces) + ". The inputs of the chosen configuration for each workflow and reference "
        "build are written next to the workflow as `example_inputs.json`, or `example_inputs.<build>.json` "
        "when it ran on more than one build. This table is for review and is not part of the repository.",
        "",
        "Choice rule, per workflow and build: the configuration with the most successful workflow runs, then "
        "the latest successful run; configurations with no successful run are never chosen. A step often ran "
        "under several configurations (for example reruns of failed samples with more resources), so check "
        "the alternatives. Override with `--choose Workflow=namespace/name`."
        + (f" Configurations whose version matches `{args.exclude_version}` were excluded." if args.exclude_version
           else "")
        + (f" Configurations for {', '.join(args.exclude_build)} were excluded." if args.exclude_build else ""),
        "",
        "| Workflow | Build | Chosen configuration | Branch run (repo copy) | Succeeded runs | Last success "
        "| Checks | Alternatives |",
        "|---|---|---|---|---:|---|---|---|",
    ]
    for w in workflows:
        groups = by_wf.get(w["name"], {})
        for old in glob.glob(os.path.join(w["dir"], "terra_config*.json")) + \
                glob.glob(os.path.join(w["dir"], "example_inputs*.json")):
            os.remove(old)
        if not groups:
            lines.append(f"| {w['name']} | | — (no configuration in the exported workspaces) | ({w['branch']}) "
                         "| | | | |")
            continue
        entry = os.path.join(w["dir"], f"{w['name']}.wdl")
        declared = womtool_inputs(entry, args.womtool)
        # A build gets its own configuration only if it ran at scale: its best configuration has at least 10%
        # of the successful runs of the workflow's best configuration. Smaller groups are usually pilots.
        best = {b: max(c["runs"].get("workflows_succeeded", 0) for c in cs) for b, cs in groups.items()}
        top = max(best.values())
        written = [b for b in sorted(groups) if best[b] > 0 and best[b] >= 0.1 * top]
        if len(written) > 1 and "unknown" in written:
            written.remove("unknown")  # a known build exists; build-less configurations are usually tests/variants
        for build in sorted(groups):
            cands = sorted(groups[build], key=rank_key, reverse=True)
            if build not in written:
                why = ("no successful run" if best[build] == 0 else
                       "build unknown, not written" if build == "unknown" else "small run, not written")
                lines.append(f"| {w['name']} | {build} | — ({why}) | | {best[build]} | | | "
                             + ", ".join(f"`{c['config']}`" for c in cands) + " |")
                continue
            chosen = next((c for c in cands if c["config"] in overrides), cands[0])
            checks = []
            version = chosen.get("method", {}).get("methodVersion")
            if version and w["branch"] and version != w["branch"]:
                checks.append("branch differs from repo copy")
            if chosen.get("current_differs_from_saved"):
                checks.append("current config edited after run")
            if chosen.get("source", "").startswith("current"):
                checks.append("no saved run config")
            if (chosen.get("run_source") or {}).get("workflow"):
                checks.append(compare_run_source(chosen, entry))
            if declared is not None:
                extra = sorted(k for k in chosen.get("inputs", {}) if k not in declared)
                if extra:
                    checks.append(f"{len(extra)} input(s) not in repo WDL: " + ", ".join(f"`{e}`" for e in extra[:3]))
            name = "example_inputs.json" if len(written) == 1 else f"example_inputs.{build}.json"
            with open(os.path.join(w["dir"], name), "w") as f:
                json.dump(example_inputs(chosen), f, indent=2)
                f.write("\n")
            runs = chosen.get("runs", {})
            alts = ", ".join(f"`{c['config']}` ({c['runs'].get('workflows_succeeded', 0)})" for c in cands
                             if c is not chosen)
            lines.append(f"| {w['name']} | {build} | `{chosen['config']}` | {version or '?'} "
                         f"({w['branch']}) | {runs.get('workflows_succeeded', 0)} | "
                         f"{runs.get('last_successful_run', '')} | {'; '.join(checks) or 'ok'} | {alts} |")

    # Both reports are for internal review (they list every configuration tried and people's Methods
    # Repository namespaces), so they go next to the git-ignored export rather than into the repository.
    export_parent = os.path.dirname(os.path.abspath(_export_dirs(args.export_dirs)[0]))
    unmatched_report = args.unmatched_report or os.path.join(export_parent, "unmatched_workflows.md")
    ulines = []
    if unmatched:
        agg: Dict[str, dict] = {}
        for r in unmatched:
            m = r.get("method", {})
            key = f"{m.get('methodPath') or (m.get('methodNamespace', '?') + '/' + m.get('methodName', '?'))}"
            a = agg.setdefault(key, {"configs": 0, "ok": 0, "last": "", "versions": set(), "source": False,
                                     "repo": m.get("sourceRepo", "?")})
            a["configs"] += 1
            a["ok"] += r.get("runs", {}).get("workflows_succeeded", 0)
            a["last"] = max(a["last"], r.get("runs", {}).get("last_successful_run", ""))
            a["versions"].add(str(m.get("methodVersion", "?")))
            a["source"] |= bool(r.get("method_source") or (r.get("run_source") or {}).get("workflow"))
        ulines += ["# Workflows run in the exported workspaces but not in this repository", "",
                  "Decide which of these belong to Phase 1. For each, the export holds the WDL that ran "
                  "when `WDL fetched` is yes (`method_source` or `run_source` in its configuration files).", "",
                  "| Method | Source | Versions run | Configurations | Succeeded runs | Last success | WDL fetched |",
                  "|---|---|---|---:|---:|---|---|"]
        for key, a in sorted(agg.items(), key=lambda kv: -kv[1]["ok"]):
            versions = ", ".join(sorted(a["versions"])[:3]) + (" …" if len(a["versions"]) > 3 else "")
            ulines.append(f"| {key} | {a['repo']} | {versions} | {a['configs']} | {a['ok']} | {a['last']} | "
                          f"{'yes' if a['source'] else 'no'} |")
        with open(unmatched_report, "w") as f:
            f.write("\n".join(ulines) + "\n")

    report = args.report or os.path.join(export_parent, "chosen_configurations.md")
    with open(report, "w") as f:
        f.write("\n".join(lines) + "\n")
    n = sum(len(v) for g in by_wf.values() for v in g.values())
    print(f"Matched {n} configurations to {len(by_wf)} workflows; {len(unmatched)} unmatched; "
          f"{excluded} excluded. Review {report}" + (f" and {unmatched_report}." if unmatched else "."))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("export", help="run on a Terra VM")
    e.add_argument("--workspace", action="append", required=True, metavar="NAMESPACE/NAME",
                   help="workspace to export; repeat for several (required)")
    e.add_argument("--out", default="terra_config_export",
                   help="output directory; one subdirectory per workspace")
    e.add_argument("--api", default=DEFAULT_API)
    e.add_argument("--orch-api", default=DEFAULT_ORCH_API)
    e.add_argument("--no-source", dest="fetch_source", action="store_false",
                   help="skip fetching the WDL that ran (faster)")

    o = sub.add_parser("organize", help="run locally from the repository root")
    o.add_argument("export_dirs", nargs="+")
    o.add_argument("--wdl-root", default="phase1/terra/wdl")
    o.add_argument("--choose", action="append", default=[], metavar="WORKFLOW=NAMESPACE/NAME",
                   help="override the chosen configuration for a workflow and build")
    o.add_argument("--womtool", help="womtool binary or jar, to check input names against the repo WDLs")
    o.add_argument("--public-bucket", action="append", default=[],
                   help="treat this bucket as public when re-checking the export")
    o.add_argument("--report", metavar="PATH",
                   help="where to write the review table of choices (default: next to the export)")
    o.add_argument("--unmatched-report", metavar="PATH",
                   help="where to list workflows not in the repository (default: next to the export)")
    o.add_argument("--exclude-build", action="append", default=[], choices=["GRCh38", "CHM13", "unknown"],
                   help="ignore configurations for this reference build; repeat for several")
    o.add_argument("--exclude-version", metavar="REGEX",
                   help="ignore configurations whose branch or snapshot matches, e.g. 'sl_v2' for Phase 2")

    args = p.parse_args(argv)
    return cmd_export(args) if args.cmd == "export" else cmd_organize(args)


if __name__ == "__main__":
    sys.exit(main())
