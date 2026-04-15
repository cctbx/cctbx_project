"""
Failure Diagnoser — prompt construction, markdown sanitiser, and HTML report
builder for the diagnosable-terminal error feature.

Kept separate from error_analyzer.py so the prompt template and HTML template
can be iterated and tested in isolation without touching detection logic.

Public API:
    build_diagnosis_prompt(error_type, error_text, program, log_tail) -> str
    _strip_llm_markdown(text) -> str
    build_diagnosis_html(description, error_excerpt, diagnosis_text,
                         program, cycle) -> str
"""

from __future__ import absolute_import, division, print_function

import html as _html_module
import re


# =============================================================================
# LLM PROMPT
# =============================================================================

def build_diagnosis_prompt(error_type, error_text, program, log_tail):
    """
    Build the LLM prompt for diagnosing a terminal program failure.

    The prompt requests plain-text output structured as three labelled
    paragraphs so the response is useful both in the HTML report and when
    read directly from the CLI log.

    The YAML diagnosis_hint is injected as a crystallographer's prior,
    steering the LLM toward the correct domain without hardcoding the answer.

    Args:
        error_type:  Key from diagnosable_errors.yaml
                     (e.g. 'crystal_symmetry_mismatch').
        error_text:  The specific error lines extracted from the failing
                     program's output.
        program:     The PHENIX program that failed (e.g. 'phenix.refine').
        log_tail:    Last section of the failing program's log file
                     (already size-capped by the caller).

    Returns:
        str: Prompt string ready for call_llm_simple().
    """
    # Import here to avoid circular imports (failure_diagnoser -> error_analyzer
    # -> yaml; error_analyzer never imports failure_diagnoser).
    # Try the production path first (libtbx-installed PHENIX), then fall back
    # to a direct relative import so tests work without a full PHENIX install.
    try:
        from libtbx.langchain.agent.error_analyzer import get_diagnosis_detector
    except ImportError:
        import os as _os
        import sys as _sys
        _agent_dir = _os.path.dirname(_os.path.abspath(__file__))
        if _agent_dir not in _sys.path:
            _sys.path.insert(0, _agent_dir)
        from error_analyzer import get_diagnosis_detector
    hint = get_diagnosis_detector().get_hint(error_type)

    # Truncate log_tail here as a second safety net (first cap is in ai_agent.py
    # before encoding; this guards against any future direct calls).
    log_section = (log_tail[-3000:] if log_tail else "(not available)")

    return (
        "You are an expert crystallographer helping a user understand why\n"
        "their PHENIX structure determination job failed.\n"
        "\n"
        "FAILED PROGRAM: {program}\n"
        "ERROR TYPE: {error_type}\n"
        "\n"
        "ERROR MESSAGE:\n"
        "{error_text}\n"
        "\n"
        "PROGRAM LOG (last section):\n"
        "{log_section}\n"
        "\n"
        "CRYSTALLOGRAPHER'S NOTES ON THIS ERROR TYPE:\n"
        "{hint}\n"
        "\n"
        "Please provide a concise, helpful diagnosis. Use plain text only --\n"
        "no markdown headers, no bullet points, no asterisks. Structure your\n"
        "response as three clearly labelled paragraphs:\n"
        "\n"
        "WHAT WENT WRONG\n"
        "[One or two sentences explaining the error in plain terms.]\n"
        "\n"
        "MOST LIKELY CAUSE\n"
        "[The most probable root cause given the error message and log.]\n"
        "\n"
        "HOW TO FIX IT\n"
        "[Specific, actionable steps. Name specific PHENIX programs or\n"
        "parameters where appropriate. Give the most likely fix first.]\n"
        "\n"
        "Keep the total response under 200 words. Do not repeat the raw\n"
        "error message verbatim."
    ).format(
        program=program or "unknown",
        error_type=error_type or "unknown",
        error_text=error_text or "(not available)",
        log_section=log_section,
        hint=hint if hint else "(none)",
    )


# =============================================================================
# MARKDOWN SANITISER
# =============================================================================

def _strip_llm_markdown(text):
    """
    Remove markdown formatting characters from LLM output.

    Applied server-side immediately after call_llm_simple() returns, before
    the text is transmitted back to the client.  The same clean plain text
    goes into both the HTML report (via white-space: pre-wrap) and the CLI log.

    Why server-side?  The client receives a single clean string that works
    everywhere; no per-surface rendering decisions are needed.

    Handles:
        ## Header / # Header   -> remove leading #, keep text
        **bold** / *italic*    -> remove asterisks, keep inner text
        __underline__          -> remove underscores (rare but seen)

    Args:
        text: Raw LLM response string.

    Returns:
        Cleaned string, or the original if text is falsy.
    """
    if not text:
        return text

    # 1. Strip leading # header markers (e.g. "## WHAT WENT WRONG")
    text = re.sub(r'^\s*#{1,6}\s*', '', text, flags=re.MULTILINE)

    # 2. Remove bold/italic markers, keeping the inner text
    #    **bold** or *italic*  →  the inner text
    text = re.sub(r'\*{1,2}([^*\n]+)\*{1,2}', r'\1', text)

    # 3. Remove underline/italic markers, keeping inner text (rare but seen).
    #    __text__ or _text_  →  the inner text.
    #
    #    CRITICAL: Only match when the opening _ is NOT preceded by a word
    #    character and the closing _ is NOT followed by a word character.
    #    This prevents eating underscores inside PHENIX parameter names
    #    (e.g. set_unit_cell_from_map=True) or filenames (nsf-d2_noligand.pdb).
    #    Real markdown italic always has non-word context around the delimiters.
    text = re.sub(r'(?<!\w)_{1,2}([^_\n]+)_{1,2}(?!\w)', r'\1', text)

    # 4. Collapse any runs of 3+ blank lines introduced by the above
    text = re.sub(r'\n{3,}', '\n\n', text)

    return text.strip()


# =============================================================================
# HTML REPORT BUILDER
# =============================================================================

def build_diagnosis_html(description, error_excerpt, diagnosis_text,
                          program, cycle,
                          html_path=None, job_name=None, working_dir=None):
    """
    Build a self-contained HTML diagnosis report.

    The diagnosis_text is expected to be plain text (already stripped of
    markdown by _strip_llm_markdown on the server).  It is placed inside a
    white-space: pre-wrap div so the three-paragraph structure is preserved
    without any markdown-to-HTML conversion.

    All user-supplied strings are passed through html.escape() before
    insertion so the report is safe even if the error excerpt or LLM output
    contains angle brackets or ampersands.

    The template has no external dependencies (no CDN, no external CSS/JS)
    so it renders correctly when opened directly from a local file path.

    Args:
        description:    Human-readable error type label from the YAML.
        error_excerpt:  The specific error lines from the failing program.
        diagnosis_text: Clean plain-text diagnosis from the LLM or fallback.
        program:        The PHENIX program that failed.
        cycle:          The agent cycle number at which the failure occurred.
        html_path:      Full path where this HTML file will be saved.
        job_name:       Human-readable job/run name (e.g. basename of working dir).
        working_dir:    Working directory of the ai_agent run.

    Returns:
        str: Complete, self-contained HTML document.
    """
    # Escape all user-supplied content before embedding in HTML.
    safe_desc    = _html_module.escape(str(description or "Unknown error"))
    safe_excerpt = _html_module.escape(str(error_excerpt or ""))
    safe_program = _html_module.escape(str(program or "unknown"))
    safe_diag    = _html_module.escape(str(diagnosis_text or
                                          "No diagnosis available."))
    safe_cycle   = _html_module.escape(str(cycle))
    safe_path    = _html_module.escape(str(html_path or ""))
    safe_job     = _html_module.escape(str(job_name or ""))
    safe_wdir    = _html_module.escape(str(working_dir or ""))

    # Build optional job-context row for the meta bar
    job_meta_html = ""
    if safe_job:
        job_meta_html += (
            "\n    &nbsp;|&nbsp;\n"
            "    Job: <strong>{safe_job}</strong>"
        ).format(safe_job=safe_job)
    if safe_wdir:
        job_meta_html += (
            "\n    &nbsp;|&nbsp;\n"
            "    Working directory: <strong>{safe_wdir}</strong>"
        ).format(safe_wdir=safe_wdir)

    # Build optional "saved to" footer line
    saved_line = ""
    if safe_path:
        saved_line = (
            "\n    Saved to: <code>{safe_path}</code><br>"
        ).format(safe_path=safe_path)

    return """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>PHENIX AI: Error Diagnosis</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica,
                   Arial, sans-serif;
      max-width: 820px;
      margin: 32px auto;
      padding: 0 20px 40px;
      line-height: 1.65;
      color: #222;
      background: #fff;
    }}
    h1 {{
      color: #c0392b;
      font-size: 1.35em;
      margin-bottom: 4px;
      border-bottom: 2px solid #c0392b;
      padding-bottom: 6px;
    }}
    h2 {{
      font-size: 1.0em;
      color: #444;
      margin-top: 1.6em;
      margin-bottom: 0.4em;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }}
    .meta {{
      color: #666;
      font-size: 0.88em;
      margin-bottom: 1.4em;
    }}
    .error-box {{
      background: #fdf3f2;
      border-left: 4px solid #c0392b;
      padding: 10px 14px;
      font-family: 'SF Mono', 'Consolas', 'Menlo', monospace;
      font-size: 0.85em;
      white-space: pre-wrap;
      word-break: break-word;
      border-radius: 0 4px 4px 0;
    }}
    .diagnosis {{
      background: #f4f8fd;
      border-left: 4px solid #2980b9;
      padding: 14px 18px;
      white-space: pre-wrap;
      font-size: 0.95em;
      border-radius: 0 4px 4px 0;
      line-height: 1.7;
    }}
    hr {{
      border: none;
      border-top: 1px solid #e0e0e0;
      margin: 1.8em 0;
    }}
    .footer {{
      color: #999;
      font-size: 0.82em;
      margin-top: 2em;
    }}
    code {{
      font-family: 'SF Mono', 'Consolas', 'Menlo', monospace;
      font-size: 0.92em;
      word-break: break-all;
    }}
  </style>
</head>
<body>
  <h1>&#9888; PHENIX AI Agent &mdash; Error diagnosis</h1>
  <div class="meta">
    Program: <strong>{safe_program}</strong>
    &nbsp;|&nbsp;
    Cycle: <strong>{safe_cycle}</strong>{job_meta_html}
  </div>

  <h2>Error</h2>
  <p><strong>{safe_desc}</strong></p>
  <div class="error-box">{safe_excerpt}</div>

  <hr>

  <h2>AI Diagnosis</h2>
  <div class="diagnosis">{safe_diag}</div>

  <hr>
  <p class="footer">
    Generated by the PHENIX AI Agent.{saved_line}
  </p>
</body>
</html>""".format(
        safe_program=safe_program,
        safe_cycle=safe_cycle,
        safe_desc=safe_desc,
        safe_excerpt=safe_excerpt,
        safe_diag=safe_diag,
        job_meta_html=job_meta_html,
        saved_line=saved_line,
    )
