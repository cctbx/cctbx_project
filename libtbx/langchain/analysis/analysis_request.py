"""Build the single analysis request that replaces the pipeline.

Step 2 of IMPLEMENTATION.md. The pipeline this replaces did:

    truncate -> one chunk -> ~1 KB summary -> retrieve -> analyse

Five measured defects lived in those arrows:

  1  the analysis step read a ~1 KB summary, not the log
     (6 of 7 logs retained under 10%)
  2  `_custom_log_chunker` dropped everything after `Citations`
     (7 of 20 logs; 59% of a 304 KB log)
  4  a fabricated program name went into the prompt
     ("Log file for phenix.xtriage_1", 20 of 20)
  5  retrieval ran against a constant string, because the query was
     built from a scrape that always returned '' and a field that was
     always ''
  6  the summariser emitted MORE than it was given, reproducibly:
     382,721 chars from 156,951 (one line 4,328 times, appearing twice
     in the source), and 2,026,595 from 312,011 (a single line of
     1.88 million characters)

**None of them can occur here, because none of those components exists
in this path.** That is the argument for the replacement -- not that a
simpler design is nicer, but that four of five defects are removed by
construction and the fifth (misidentification) is fixed in
`program_identity`.

WHAT THIS MODULE GUARANTEES

* the log is embedded **whole, once, byte for byte**, and a SHA-256 of
  it travels with the request so a caller can verify before sending;
* the header names a program only when one was resolved, and says so
  plainly when none was;
* the payload cannot be larger than the log plus a bounded instruction
  block -- defect 6 turned inside out into an assertion.

WHAT IT DOES NOT DO

No LLM call. This builds the request; the caller sends it. That split
is deliberate: every guarantee above is testable without a key, a
network or a provider, which is why `tst_analysis_request` runs in the
ordinary suite.
"""
from __future__ import absolute_import, division, print_function

import hashlib

#: Delimiters. Explicit rather than implicit so the log block can be
#: located and hashed by anything downstream, including a test.
LOG_OPEN = "<<<PHENIX_LOG_BEGIN>>>"
LOG_CLOSE = "<<<PHENIX_LOG_END>>>"

INSTRUCTIONS = """\
You are an expert crystallographer reading a Phenix program log.

Report, in this order:
  1. What the run did, and on what inputs.
  2. Whether it worked, with the metrics the log states.
  3. Anything that went wrong: warnings, errors, unexpected values.
  4. What to do next, naming the Phenix programs and their inputs.

Rules:
  * Use only what the log states. Do not supply values it does not give.
  * If the log does not say something, say that it does not say it.
  * Be concise. Prefer the log's own numbers to prose about them.
"""

#: A candidate intervention, deliberately NOT in the default prompt.
#:
#: The first version of INSTRUCTIONS contained:
#:
#:     "do not recommend molecular replacement as the primary route on
#:      a dataset with a usable anomalous signal"
#:
#: and Step 2's acceptance test is "must not recommend molecular
#: replacement as the primary route". **That is teaching to the test.**
#: With the rule in the prompt, passing proves the model can follow an
#: instruction, not that the whole log reached it and was reasoned
#: about -- which is the only thing Step 2 needs to establish.
#:
#: The rule may well improve real reports. If so it can be measured as
#: an intervention, on a baseline established without it, and never
#: with the same statement serving as both treatment and test.
PHASING_HINT = """\
  * Where the data support experimental phasing, say so; do not
    recommend molecular replacement as the primary route on a dataset
    with a usable anomalous signal.
"""

NO_PROGRAM = ("The program that produced this log is not stated in the "
              "log and could not be determined. Do not guess it.")


class AnalysisRequest(object):
    """What the caller sends, plus what a test needs to verify it."""

    __slots__ = ("payload", "header", "log_sha256", "log_chars",
                 "payload_chars", "program", "program_source")

    def __init__(self, payload, header, log_sha256, log_chars,
                 program, program_source):
        self.payload = payload
        self.header = header
        self.log_sha256 = log_sha256
        self.log_chars = log_chars
        self.payload_chars = len(payload)
        self.program = program
        self.program_source = program_source

    @property
    def overhead(self):
        """Characters added to the log. Bounded by construction; the
        component this replaces added 1.7 million on one input."""
        return self.payload_chars - self.log_chars

    def __repr__(self):
        return ("AnalysisRequest(program=%r, source=%r, log=%d, payload=%d)"
                % (self.program, self.program_source, self.log_chars,
                   self.payload_chars))


def sha256(text):
    return hashlib.sha256(text.encode("utf-8", "replace")).hexdigest()


def build_analysis_request(log_text, identity=None, instructions=None):
    """Build the request for one long-context analysis call.

    `identity` is a `program_identity.Identity` or None. When it carries
    no name, **no header is written** -- the shipped code wrote
    "Log file for phenix.<derived>" unconditionally and named a
    nonexistent program on 20 of 20 corpus logs.
    """
    log_text = log_text or ""
    name = getattr(identity, "name", None)
    source = getattr(identity, "source", "none")

    if name:
        header = "This is a log file from %s." % name
        if getattr(identity, "authoritative", False) is False:
            header += (" (inferred from the file name; the log does not "
                       "state it)")
    else:
        header = ""

    parts = [instructions or INSTRUCTIONS, ""]
    parts.append(header if header else NO_PROGRAM)
    parts.append("")
    parts.append(LOG_OPEN)
    parts.append(log_text)
    parts.append(LOG_CLOSE)

    payload = "\n".join(parts)
    return AnalysisRequest(payload, header, sha256(log_text), len(log_text),
                           name, source)


def extract_log_block(payload):
    """Recover the embedded log. Lets a caller verify the hash on the
    thing it is actually about to send, rather than on what it believes
    it embedded."""
    start = payload.find(LOG_OPEN)
    end = payload.find(LOG_CLOSE)
    if start == -1 or end == -1:
        return ""
    return payload[start + len(LOG_OPEN) + 1:end - 1]


def verify(request):
    """Check the request before it goes out.

    Returns (ok, [problems]). Cheap, and the only guard between a
    construction bug and a user-facing report built on the wrong text.
    """
    problems = []
    block = extract_log_block(request.payload)
    if sha256(block) != request.log_sha256:
        problems.append("the embedded log does not hash to the source")
    if request.payload.count(LOG_OPEN) != 1:
        problems.append("the log block is not present exactly once")
    if request.overhead < 0:
        problems.append("payload is smaller than the log")
    if request.program is None and request.program_source != "none":
        problems.append("no program but a source is claimed")
    return (not problems), problems
