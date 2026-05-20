"""
Sandbox test: verify the Python dependency environment for the
PHENIX AI agent.

Purpose
-------
Run this as a pre-flight check to confirm the current Python env
has everything the ai_agent needs.  Unlike the K_A..K_G tests
which verify code behavior, this test verifies the environment:
whether each required package can be imported, and (for packages
known to have transitive version conflicts) whether they import
cleanly without raising during init.

Categories
----------
1. REQUIRED (must succeed for any use of ai_agent):
   - langchain_core
   - langchain_community
   - yaml (pyyaml)

2. PROVIDER (at least one must succeed; user picks via config):
   - google: langchain_google_genai + EITHER google.genai (new)
     OR google.generativeai (deprecated; will be removed by Google)
   - openai: langchain_openai (+ openai)
   - anthropic: langchain_anthropic (+ anthropic)
   - ollama: langchain_ollama

3. RAG-OPTIONAL (test SKIPS, doesn't FAIL, if missing — Section G
   makes the code gracefully degrade when these are absent):
   - chromadb
   - langchain_chroma
   - langchain_cohere + cohere

4. DOC-LOADING (only needed to rebuild the RAG database; users
   running pre-built db don't need these):
   - pypdf
   - unstructured + bs4
   - langchain_text_splitters
   - tqdm

5. RENDERING (cosmetic only):
   - markdown_it
   - linkify_it

6. RUNTIME PROBES (catch "installed but broken" failure modes):
   - Chroma: `from langchain_chroma import Chroma` — catches the
     v118.G protobuf-version conflict.  If both chromadb and
     langchain_chroma are pip-installed but this probe fails,
     suggest the opentelemetry-proto upgrade documented in
     v118.6.1's install hint.
   - FlashrankRerank: needed at runtime by analyze_log_summary()
     and query_docs() via create_reranking_retriever().  If
     chromadb is available but flashrank isn't, RAG queries will
     crash — surface this as FAIL so the user installs flashrank
     before first use.

Test outcomes
-------------
- PASS:  all REQUIRED + at least one PROVIDER + (RAG probe OK if
         chromadb/langchain_chroma installed)
- FAIL:  any REQUIRED missing, OR all PROVIDERs missing, OR RAG
         probe fails (chromadb+langchain_chroma installed but
         can't load Chroma — actionable env issue)
- WARN:  any non-REQUIRED package missing (informational)

Run
---
    libtbx.python tst_dependencies.py
    # or
    libtbx.python tst_dependencies.py --strict   # also fail on warnings
"""

from __future__ import absolute_import, division, print_function

import importlib
import sys


# ---------------------------------------------------------------
# Package classifications
# ---------------------------------------------------------------

# Format: list of (pip_name, import_name, description) tuples.
# Some packages have different pip and import names (e.g.
# `pyyaml` → `yaml`, `beautifulsoup4` → `bs4`).

REQUIRED = [
    ("langchain-core",      "langchain_core",      "Core LangChain runtime"),
    ("langchain-community", "langchain_community", "Document loaders, retrievers"),
    ("pyyaml",              "yaml",                "knowledge/programs.yaml loading"),
]

PROVIDERS = [
    # (pip_name, primary_import_name, description, [(extra_pip, extra_import), ...])
    # `extra_import` is given EXPLICITLY (not derived) because some packages have
    # dotted import paths or names that don't match pip name + s/-/_/.
    #
    # GOOGLE: google.generativeai is being deprecated in favor of google.genai
    # (https://github.com/google-gemini/deprecated-generative-ai-python).  Code
    # in agent/directive_extractor.py already tries the new package first and
    # falls back.  The env check treats EITHER package as satisfying the extra:
    # see _probe_google_extras() below.
    ("langchain-google-genai", "langchain_google_genai", "Google Gemini provider",
        "GOOGLE_SPECIAL"),
    ("langchain-openai",       "langchain_openai",       "OpenAI provider",
        [("openai", "openai")]),
    ("langchain-anthropic",    "langchain_anthropic",    "Anthropic provider",
        [("anthropic", "anthropic")]),
    ("langchain-ollama",       "langchain_ollama",       "Local Ollama provider",
        []),
]

RAG_OPTIONAL = [
    ("langchain-chroma", "langchain_chroma", "RAG vector store wrapper"),
    ("chromadb",         "chromadb",         "RAG vector DB engine"),
    ("langchain-cohere", "langchain_cohere", "Cohere reranker (alternative)"),
    ("cohere",           "cohere",           "Cohere API client"),
    # flashrank moved to runtime probe — see _probe_flashrank_runtime() and
    # _report_flashrank_runtime_probe().  It's needed at runtime by both
    # analyze_log_summary() and query_docs() via create_reranking_retriever(),
    # so missing-flashrank surfaces as a runtime failure, not a warning.
]

DOC_LOADING_OPTIONAL = [
    ("pypdf",                    "pypdf",                    "PDF document loader"),
    ("unstructured",             "unstructured",             "HTML document loader"),
    ("beautifulsoup4",           "bs4",                      "HTML parser (unstructured backend)"),
    ("langchain-text-splitters", "langchain_text_splitters", "Document chunking"),
    ("tqdm",                     "tqdm",                     "Progress bars"),
]

RENDERING_OPTIONAL = [
    ("markdown-it-py", "markdown_it", "Markdown rendering"),
    ("linkify-it-py",  "linkify_it",  "Link detection in markdown"),
]


# ---------------------------------------------------------------
# Probe helpers
# ---------------------------------------------------------------

def _probe(import_name):
    """Try to import. Returns (ok, error_string)."""
    try:
        importlib.import_module(import_name)
        return (True, None)
    except Exception as e:
        # Catch Exception, not just ImportError: version conflicts
        # surface as TypeError, RuntimeError, etc.  This is the
        # same lesson Section G encoded — see
        # docs/DEVELOPER_GUIDE.md "Optional dependency handling".
        return (False, "%s: %s" % (type(e).__name__, e))


def _probe_chroma_class():
    """Special probe: actually load the Chroma class.

    chromadb's transitive opentelemetry-proto dependency can be
    pinned at an ancient version that ships _pb2.py files
    incompatible with newer protobuf.  Plain `import chromadb`
    may succeed while `from langchain_chroma import Chroma`
    raises TypeError during langchain_chroma's __init__.

    This probe catches the v118.G failure mode at env-check
    time, before any user code hits the issue.
    """
    try:
        from langchain_chroma import Chroma
        # Confirm Chroma is the expected class
        if Chroma is None:
            return (False, "langchain_chroma.Chroma is None")
        return (True, None)
    except Exception as e:
        return (False, "%s: %s" % (type(e).__name__, e))


def _probe_google_extras():
    """Probe for the Google Gemini SDK with deprecation-awareness.

    The `google-generativeai` package (Python import: `google.generativeai`)
    has been deprecated by Google in favor of `google-genai` (Python
    import: `google.genai`).  Both packages can satisfy ai_agent's
    Google provider — `directive_extractor.py` tries the new one
    first and falls back to the old one.

    Returns:
        (ok, found_pip_name, error_string)
        - ok=True if EITHER package is importable
        - found_pip_name names which one was found (or None on failure)
        - error_string describes the failure mode if neither importable
    """
    # Try the new package first
    ok_new, err_new = _probe("google.genai")
    if ok_new:
        return (True, "google-genai", None)
    # Fall back to the deprecated one
    ok_old, err_old = _probe("google.generativeai")
    if ok_old:
        return (True, "google-generativeai (deprecated; switch to google-genai)", None)
    # Neither available
    return (False, None,
            "neither google.genai (%s) nor google.generativeai (%s) is importable"
            % (err_new, err_old))


def _probe_flashrank_runtime():
    """Probe whether FlashrankRerank can actually be instantiated.

    flashrank is needed at runtime by analyze_log_summary() and
    query_docs() via create_reranking_retriever() in rag/retriever.py.
    Plain `import flashrank` may succeed while instantiation fails
    (e.g. missing native binaries, model download issues).

    Detects the analogous failure mode to v118.G's chromadb probe:
    "package imports but symbol cannot be loaded."
    """
    try:
        from langchain_community.document_compressors import FlashrankRerank
        if FlashrankRerank is None:
            return (False, "FlashrankRerank is None after import")
        # Don't instantiate — instantiation tries to download model weights.
        # Just confirm the class is constructible (has __init__).
        if not callable(FlashrankRerank):
            return (False, "FlashrankRerank is not callable")
        return (True, None)
    except Exception as e:
        return (False, "%s: %s" % (type(e).__name__, e))


# ---------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------

def _report_flashrank_runtime_probe():
    """Probe whether FlashrankRerank can be loaded for reranking.

    Returns True if usable OR if the user is clearly opting out of
    RAG entirely (chromadb missing → no RAG path will be taken).
    """
    print()
    print("=== Reranker runtime probe ===")
    # Skip if RAG itself isn't available — no point probing the reranker
    chroma_ok, _ = _probe("chromadb")
    if not chroma_ok:
        print("  [SKIP]    Reranker probe — chromadb not installed; "
              "the reranking path is unreachable in this env")
        return True
    ok, err = _probe_flashrank_runtime()
    if ok:
        print("  [ OK ]    FlashrankRerank loaded; reranking retriever usable")
        return True
    print("  [FAIL]    FlashrankRerank could not be loaded: %s" % err)
    print()
    print("  flashrank is needed at runtime by both `analyze_log_summary()`")
    print("  and `query_docs()` via `create_reranking_retriever()`.  Without")
    print("  it, every RAG query will fail.  Install:")
    print()
    print("    pip install flashrank")
    print()
    print("  Or, if you want to use RAG without reranking, code changes")
    print("  in analyzer.py/query.py would be needed to use a non-reranking")
    print("  retriever path.")
    return False


def _report_category(name, entries, *, severity):
    """Probe every package in a category and print PASS/FAIL/WARN
    for each.  Returns (passed_pip_names, failed_pip_names).
    """
    print()
    print("=== %s ===" % name)
    passed = []
    failed = []
    for entry in entries:
        if len(entry) == 4:
            pip_name, import_name, desc, _extras = entry
        else:
            pip_name, import_name, desc = entry
        ok, err = _probe(import_name)
        if ok:
            print("  [ OK ]    %s (%s) — %s" % (pip_name, import_name, desc))
            passed.append(pip_name)
        else:
            tag = "[FAIL]" if severity == "required" else \
                  "[WARN]" if severity in ("rag", "doc", "render") else \
                  "[----]"  # provider — fail-once handled by caller
            print("  %s    %s (%s) — %s" % (tag, pip_name, import_name, desc))
            print("              error: %s" % err)
            print("              install: pip install %s" % pip_name)
            failed.append(pip_name)
    return passed, failed



    """Probe every package in a category and print PASS/FAIL/WARN
    for each.  Returns (passed_pip_names, failed_pip_names).
    """
    print()
    print("=== %s ===" % name)
    passed = []
    failed = []
    for entry in entries:
        if len(entry) == 4:
            pip_name, import_name, desc, _extras = entry
        else:
            pip_name, import_name, desc = entry
        ok, err = _probe(import_name)
        if ok:
            print("  [ OK ]    %s (%s) — %s" % (pip_name, import_name, desc))
            passed.append(pip_name)
        else:
            tag = "[FAIL]" if severity == "required" else \
                  "[WARN]" if severity in ("rag", "doc", "render") else \
                  "[----]"  # provider — fail-once handled by caller
            print("  %s    %s (%s) — %s" % (tag, pip_name, import_name, desc))
            print("              error: %s" % err)
            print("              install: pip install %s" % pip_name)
            failed.append(pip_name)
    return passed, failed


def _report_chroma_runtime_probe():
    """Run the Section-G failure-mode probe and report."""
    print()
    print("=== RAG runtime probe (v118.G) ===")
    ok, err = _probe_chroma_class()
    if ok:
        print("  [ OK ]    `from langchain_chroma import Chroma` succeeds")
        return True
    # Distinguish "not installed" from "installed but broken"
    chroma_ok, _ = _probe("chromadb")
    lc_ok, _ = _probe("langchain_chroma")
    if not chroma_ok or not lc_ok:
        # Either package is just missing — not a config error;
        # already covered by RAG_OPTIONAL section above.
        print("  [SKIP]    Chroma probe — chromadb or langchain_chroma not "
              "installed; this is expected for users not using RAG")
        return True
    # Packages installed but probe failed — this is the v118.G case
    print("  [FAIL]    `from langchain_chroma import Chroma` raised: %s"
          % err)
    print()
    print("  This is the v118.G failure mode.  Both langchain_chroma and")
    print("  chromadb are installed, but a transitive dep is at an")
    print("  incompatible version.  Most likely: opentelemetry-proto is")
    print("  too old and ships _pb2.py files incompatible with the")
    print("  installed protobuf.  Try:")
    print()
    print("    pip install --upgrade opentelemetry-proto "
          "opentelemetry-exporter-otlp-proto-grpc")
    print()
    print("  Or, as a no-package-changes workaround (slower runtime):")
    print()
    print("    export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python")
    print()
    return False


# ---------------------------------------------------------------
# Main test
# ---------------------------------------------------------------

def check_environment(strict=False):
    """Run env-readiness checks.  Returns True iff env is OK.

    This is the core function — call it programmatically to get a
    boolean.  Used directly by the standalone runner below and by
    run_all_tests() (which adapts the return to the test-suite
    convention of "raise on failure").

    Args:
        strict: if True, treat WARN (missing optional packages)
                as failure too.  Default False — only REQUIRED
                + provider availability + RAG runtime probe
                determine pass/fail.
    """
    print("=" * 65)
    print("PHENIX AI agent — environment readiness check")
    print("=" * 65)
    print()
    print("Python: %s" % sys.version.replace("\n", " "))
    print()

    # 1. REQUIRED
    req_passed, req_failed = _report_category(
        "REQUIRED packages", REQUIRED, severity="required")

    # 2. PROVIDERS — at least one must work
    print()
    print("=== PROVIDER packages (need ≥1 to make LLM calls) ===")
    providers_ok = []
    for pip_name, import_name, desc, extras in PROVIDERS:
        ok, err = _probe(import_name)
        # Probe each extra by its EXPLICIT import name (do NOT derive
        # via s/-/_/ — `google-generativeai` imports as `google.generativeai`).
        #
        # GOOGLE_SPECIAL: the Google provider needs either google-genai
        # (current) or google-generativeai (deprecated).  Either suffices.
        extras_ok = True
        extras_err = None
        extras_found = None  # for display: what was actually found
        if extras == "GOOGLE_SPECIAL":
            xok, found_pkg, xerr = _probe_google_extras()
            if not xok:
                extras_ok = False
                extras_err = xerr
            else:
                extras_found = found_pkg
            display_extras = [("google-genai or google-generativeai", None)]
            install_extras = ["google-genai"]  # recommend the new one
        else:
            for extra_pip, extra_import in extras:
                xok, xerr = _probe(extra_import)
                if not xok:
                    extras_ok = False
                    extras_err = "%s (import %s): %s" % (extra_pip, extra_import, xerr)
                    break
            display_extras = list(extras)
            install_extras = [p for p, _ in extras]

        if ok and extras_ok:
            line = "  [ OK ]    %s — %s" % (pip_name, desc)
            if extras_found:
                # GOOGLE_SPECIAL path — show which one was found
                line += "  (+ %s)" % extras_found
            elif display_extras:
                line += "  (+ %s)" % ", ".join(p for p, _ in display_extras)
            print(line)
            providers_ok.append(pip_name)
        else:
            print("  [----]    %s — %s" % (pip_name, desc))
            if not ok:
                print("              error: %s" % err)
            if not extras_ok:
                print("              extras missing: %s" % extras_err)
            install_cmd = "pip install %s" % pip_name
            if install_extras:
                install_cmd += " " + " ".join(install_extras)
            print("              install: %s" % install_cmd)

    # 3. RAG-OPTIONAL
    _, _ = _report_category(
        "RAG packages (optional; Section G makes the agent degrade gracefully)",
        RAG_OPTIONAL, severity="rag")

    # 4. Runtime probe — catches v118.G's protobuf failure case
    chroma_probe_ok = _report_chroma_runtime_probe()

    # 4b. Reranker runtime probe — catches missing/broken flashrank
    flashrank_probe_ok = _report_flashrank_runtime_probe()

    # 5. DOC-LOADING
    _, doc_failed = _report_category(
        "DOC-LOADING packages (optional; only needed to rebuild docs DB)",
        DOC_LOADING_OPTIONAL, severity="doc")

    # 6. RENDERING
    _, render_failed = _report_category(
        "RENDERING packages (cosmetic only)",
        RENDERING_OPTIONAL, severity="render")

    # =============================================================
    # Verdict
    # =============================================================
    print()
    print("=" * 65)
    print("Summary")
    print("=" * 65)

    fatal = []
    if req_failed:
        fatal.append("REQUIRED packages missing: %s" %
                     ", ".join(req_failed))
    if not providers_ok:
        fatal.append("No PROVIDER package is usable — at least one of "
                     "(google, openai, anthropic, ollama) must be "
                     "installed and importable")
    if not chroma_probe_ok:
        fatal.append("RAG runtime probe failed — chromadb stack "
                     "installed but broken (see remediation above)")
    if not flashrank_probe_ok:
        fatal.append("Reranker runtime probe failed — chromadb is "
                     "available but flashrank is not, so RAG queries "
                     "will crash at runtime (see remediation above)")

    if fatal:
        print()
        print("RESULT: FAIL")
        print()
        for f in fatal:
            print("  - %s" % f)
        print()
        print("Optional packages missing (informational):")
        all_optional_missing = doc_failed + render_failed
        if all_optional_missing:
            for p in all_optional_missing:
                print("    %s" % p)
        else:
            print("    (none)")
        return False

    print()
    print("RESULT: PASS")
    print("  Required: %d/%d" % (len(req_passed), len(REQUIRED)))
    print("  Provider: %s" % ", ".join(providers_ok))
    print("  RAG runtime probe: OK")
    print("  Reranker runtime probe: OK")

    optional_missing = doc_failed + render_failed
    if optional_missing:
        print()
        print("Optional packages missing (does not block ai_agent):")
        for p in optional_missing:
            print("    %s" % p)
        if strict:
            print()
            print("RESULT (strict mode): FAIL — optional packages missing")
            return False

    print()
    print("Environment ready for PHENIX AI agent.")
    return True


def run_all_tests(strict=False):
    """Test-suite-compatible entry point.

    The Phenix test runner (run_all_tests.py:run_test_module) signals
    test failure by uncaught exception.  We adapt check_environment()'s
    boolean return value to that convention: raise AssertionError on
    FAIL so the registry records it correctly.

    Use check_environment() directly (or the standalone __main__ runner)
    if you want the boolean return without the exception side-effect.
    """
    ok = check_environment(strict=strict)
    if not ok:
        raise AssertionError(
            "Environment dependency check failed — see output above. "
            "Run `libtbx.python tst_dependencies.py` directly for the "
            "full diagnostic, or `libtbx.python tst_dependencies.py "
            "--strict` to additionally fail on optional missing.")
    return True


if __name__ == "__main__":
    strict = "--strict" in sys.argv
    success = check_environment(strict=strict)
    sys.exit(0 if success else 1)
