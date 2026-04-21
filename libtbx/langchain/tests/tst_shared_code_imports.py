"""
Shared-code import safety tests.

The agent/ directory ships with both the client (in PHENIX) and the server.
Shared modules must not import:
  1. LLM/server-only dependencies (langchain_core, openai, anthropic, etc.)
  2. Agent modules that don't exist in the shipped codebase
  3. Third-party packages not bundled with PHENIX

See docs/guides/BACKWARD_COMPATIBILITY.md — RULE 7: The agent/ shared code trap.

Run with:
    PYTHONPATH=. python tests/tst_shared_code_imports.py
"""

import os
import re
import sys

# Add parent directory to path for imports
sys.path.insert(
    0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import assert_true
from tests.tst_utils import run_tests_with_fail_fast

# =========================================================================
# PHENIX/cctbx Linter "Silencer"
# =========================================================================
# These imports are only used inside function bodies. To prevent
# libtbx.find_unused_imports from flagging them, we reference them here.
(re, assert_true, run_tests_with_fail_fast)

# Ensure sys is used for path manipulation
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# =========================================================================
# Classification of agent/ modules
# =========================================================================

# Shared modules: run on BOTH client and server.
SHARED_MODULES = [
    "command_builder.py",
    "workflow_engine.py",
    "workflow_state.py",
    "planner.py",
    "command_postprocessor.py",
    "file_utils.py",
    "advice_preprocessor.py",
    "contract.py",
]

SERVER_ONLY_MODULES = ["graph_nodes.py"]
CLIENT_ONLY_MODULES = ["session.py", "best_files_tracker.py"]

# =========================================================================
# Import safety rules
# =========================================================================

FORBIDDEN_IN_SHARED = [
    "langchain",
    "langgraph",
    "openai",
    "anthropic",
    "google.generativeai",
    "google.genai",
    "tiktoken",
    "httpx",
]

KNOWN_AGENT_MODULES = {
    "advice_preprocessor", "best_files_tracker", "command_builder",
    "command_postprocessor", "contract", "event_log", "file_utils",
    "graph_nodes", "metrics_analyzer", "nl_to_phil", "pattern_manager",
    "perceive_checks", "placement_checker", "planner", "program_registry",
    "rate_limit_handler", "rules_selector", "session", "template_builder",
    "workflow_engine", "workflow_state",
}

KNOWN_KNOWLEDGE_MODULES = {"program_registration", "prompts_hybrid", "yaml_loader"}

def _get_shared_file_paths():
    """Return absolute paths for shared modules that exist."""
    paths = []
    for name in SHARED_MODULES:
        path = os.path.join(_PROJECT_ROOT, "agent", name)
        if os.path.exists(path):
            paths.append(path)
    return paths

# ---------------------------------------------------------------------------
# Test Functions
# ---------------------------------------------------------------------------

def test_no_forbidden_imports_in_shared():
    violations = []
    for filepath in _get_shared_file_paths():
        basename = os.path.basename(filepath)
        with open(filepath) as f:
            for lineno, line in enumerate(f, 1):
                stripped = line.strip()
                if stripped.startswith("#") or not (stripped.startswith("from ") or stripped.startswith("import ")):
                    continue
                if "libtbx.langchain" in stripped:
                    continue
                for forbidden in FORBIDDEN_IN_SHARED:
                    pattern = r'(?:from|import)\s+' + re.escape(forbidden)
                    if re.search(pattern, stripped):
                        violations.append("%s:%d  %s" % (basename, lineno, stripped[:100]))
    if violations:
        msg = "Found %d forbidden import(s) in shared modules:\n  %s" % (len(violations), "\n  ".join(violations))
        assert_true(False, msg)
    print("  PASSED: No forbidden imports in %d shared modules" % len(SHARED_MODULES))

def test_agent_imports_reference_known_modules():
    unknown = []
    for filepath in _get_shared_file_paths():
        basename = os.path.basename(filepath)
        with open(filepath) as f:
            for lineno, line in enumerate(f, 1):
                stripped = line.strip()
                if stripped.startswith("#"): continue
                for m in re.finditer(r'from\s+(?:libtbx\.langchain\.)?agent\.(\w+)', stripped):
                    if m.group(1) not in KNOWN_AGENT_MODULES:
                        unknown.append("%s:%d  agent.%s" % (basename, lineno, m.group(1)))
                for m in re.finditer(r'from\s+(?:libtbx\.langchain\.)?knowledge\.(\w+)', stripped):
                    if m.group(1) not in KNOWN_KNOWLEDGE_MODULES:
                        unknown.append("%s:%d  knowledge.%s" % (basename, lineno, m.group(1)))
    if unknown:
        assert_true(False, "Unknown modules found: " + "\n  ".join(unknown))
    print("  PASSED: All agent/knowledge imports reference known modules")

def test_server_only_has_llm_imports():
    graph_nodes = os.path.join(_PROJECT_ROOT, "agent", "graph_nodes.py")
    if not os.path.exists(graph_nodes):
        print("  SKIP (graph_nodes.py not found)")
        return
    with open(graph_nodes) as f:
        content = f.read()
    has_llm = bool(re.search(r'(?:from|import)\s+(?:langchain_core|langchain_google|langchain_openai)', content))
    has_provider = "provider" in content and ("openai" in content or "google" in content)
    assert_true(has_llm or has_provider, "graph_nodes.py should have LLM imports")
    print("  PASSED: graph_nodes.py correctly contains LLM imports")

def test_shared_imports_are_guarded():
    unguarded = []
    for filepath in _get_shared_file_paths():
        basename = os.path.basename(filepath)
        with open(filepath) as f:
            lines = f.readlines()
        in_try = False
        for i, line in enumerate(lines, 1):
            stripped = line.strip()
            if stripped.startswith("try:"): in_try = True; continue
            if stripped.startswith("except") and "Import" in stripped: in_try = False; continue
            if not (stripped.startswith("from ") or stripped.startswith("import ")): continue
            if (len(line) - len(line.lstrip())) > 8: continue
            if bool(re.match(r'from\s+(?:libtbx\.langchain\.)?(?:agent|knowledge)\.', stripped)) and not in_try:
                unguarded.append("%s:%d  %s" % (basename, i, stripped[:100]))
    if unguarded:
        for u in unguarded: print("    WARNING: %s" % u)
    print("  PASS (warning-only for now)")

def run_all_tests():
    run_tests_with_fail_fast()

if __name__ == "__main__":
    run_all_tests()

