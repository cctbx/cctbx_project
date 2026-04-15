"""
Base class for program-specific validators.

Validators handle program-specific quirks and requirements:
- Add required parameters (e.g., R-free flags for phenix.refine)
- Fix common syntax errors
- Check prerequisites

To create a new validator:
1. Subclass ProgramValidator
2. Implement program_name property
3. Override prevalidate/postvalidate as needed
4. Place in validation/phenix_yourprogram.py (auto-discovered)

Example:
    class MyProgramValidator(ProgramValidator):
        @property
        def program_name(self) -> str:
            return "phenix.my_program"

        def prevalidate(self, command, project_state):
            # Add required parameter if missing
            if 'required_param=' not in command:
                command += ' required_param=True'
            return command
"""
from __future__ import absolute_import, division, print_function

from abc import ABC, abstractmethod
from typing import Optional
assert Optional is not None
assert abstractmethod is not None

class ProgramValidator(ABC):
    """
    Abstract base class for program-specific validation.

    Subclasses must implement:
    - program_name: Which program this validates

    Optional overrides:
    - prevalidate(): Modify command before syntax validation
    - postvalidate(): Modify command after syntax validation
    - get_common_errors(): Return known error->fix mappings
    """

    @property
    @abstractmethod
    def program_name(self) -> str:
        """
        Name of the Phenix program this validator handles.

        Must match how the program is invoked (e.g., 'phenix.refine').

        Returns:
            str: Program name (e.g., 'phenix.refine', 'phenix.phaser')
        """
        pass

    @property
    def description(self) -> str:
        """
        Human-readable description of what this validator does.

        Returns:
            str: Description for documentation
        """
        return f"Validator for {self.program_name}"

    def prevalidate(self, command: str, project_state: dict) -> str:
        """
        Modify command BEFORE syntax validation.

        Use this to:
        - Add required parameters that must be present
        - Fix known issues preemptively
        - Check prerequisites

        Args:
            command: The command to validate
            project_state: Current project state (for context)

        Returns:
            str: The (potentially modified) command

        Example:
            def prevalidate(self, command, project_state):
                # Add R-free generation if MTZ lacks it
                if 'r_free_flags.generate' not in command:
                    if not self._mtz_has_rfree(command):
                        command += ' xray_data.r_free_flags.generate=True'
                return command
        """
        return command

    def postvalidate(self, command: str, project_state: dict) -> str:
        """
        Modify command AFTER syntax validation passes.

        Use this to:
        - Add optional improvements
        - Fine-tune parameters based on project state

        Args:
            command: The validated command
            project_state: Current project state

        Returns:
            str: The (potentially modified) command
        """
        return command

    def get_common_errors(self) -> dict:
        """
        Return known error patterns and their fixes.

        Returns:
            dict: Mapping of error pattern -> fix information

        Example:
            def get_common_errors(self):
                return {
                    "Ambiguous parameter.*twin_law": {
                        "fix": "Use full path: refinement.main.twin_law=...",
                        "example": "refinement.main.twin_law=-h,-k,l"
                    }
                }
        """
        return {}

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}(program='{self.program_name}')>"
