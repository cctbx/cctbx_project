"""
Validator for phenix.refine.

Handles phenix.refine-specific validation:
- R-free flag generation
- Twin law parameter paths
"""
from __future__ import absolute_import, division, print_function

from libtbx.langchain.validation.base import ProgramValidator
from libtbx.langchain.validation.mtz_utils import add_rfree_generation_if_needed


class RefineValidator(ProgramValidator):
    """
    Validator for phenix.refine.

    Handles:
    - R-free flag generation when MTZ lacks test set
    - Twin law parameter path fixes
    """

    @property
    def program_name(self) -> str:
        return "phenix.refine"

    @property
    def description(self) -> str:
        return "Validator for phenix.refine - handles R-free flags and twin laws"

    def prevalidate(self, command: str, project_state: dict) -> str:
        """
        Add R-free generation if MTZ file lacks test set.
        """
        command = add_rfree_generation_if_needed(command, self.program_name)
        return command

    def get_common_errors(self) -> dict:
        """
        Return known error patterns and fixes for phenix.refine.
        """
        return {
            "Ambiguous parameter.*twin_law": {
                "fix": "Use full path: refinement.main.twin_law=...",
                "example": "refinement.main.twin_law=-h,-k,l"
            },
            "No array of R-free flags found": {
                "fix": "Add: xray_data.r_free_flags.generate=True",
                "always_works": True
            },
            "Ambiguous parameter.*number_of_macro_cycles": {
                "fix": "Use full path: refinement.main.number_of_macro_cycles=...",
                "example": "refinement.main.number_of_macro_cycles=5"
            }
        }
