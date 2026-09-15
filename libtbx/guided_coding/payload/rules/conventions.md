# Code and documentation conventions for this module (cctbx style)

Grounded in files on cctbx_project master (for example `libtbx/str_utils.py`, `libtbx/test_utils/__init__.py`), 2026-09-10. Confirm against the surrounding code of any file you touch; the local code wins over this file.

## Code

- Two-space indentation. Match the file you are editing exactly.
- At the top of every Python file (after any header comment): `from __future__ import absolute_import, division, print_function`.
- Plain functions and small classes in the local style; no new frameworks, no dependency additions, no broad style retrofits. Change only what the plan names.
- Follow the naming and layout of the surrounding module, not generic Python fashion. When the surrounding code and PEP 8 disagree, the surrounding code wins.

## Documentation of new code

- Every new public function, class, and method gets a docstring: one line saying what it does, then parameters and return value where they are not obvious. Match the docstring style already used in the module; where the module has none, use plain one-paragraph docstrings.
- A new module gets a short header comment stating its purpose.
- Changed behavior of an existing public function gets its docstring updated in the same change.
- The reviewer flags new public code that lacks these; missing documentation is a finding, not a nit to skip.
