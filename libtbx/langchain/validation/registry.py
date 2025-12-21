"""
Validator registry for automatic discovery and lookup.

This module maintains a registry of all program validators,
allowing easy lookup by program name.

Usage:
    from libtbx.langchain.validation import get_validator

    validator = get_validator('phenix.refine')
    if validator:
        command = validator.prevalidate(command, project_state)
"""
from __future__ import absolute_import, division, print_function

from libtbx.langchain.validation.base import ProgramValidator
from libtbx.langchain.validation.phenix_refine import RefineValidator


# =============================================================================
# Validator Registry
# =============================================================================

# Register all validators here
_VALIDATORS = {
    'phenix.refine': RefineValidator(),
    # Add more validators as they are created:
    # 'phenix.phaser': PhaserValidator(),
    # 'phenix.xtriage': XtriageValidator(),
}


def get_validator(program_name: str) -> ProgramValidator:
    """
    Get the validator for a specific program.

    Args:
        program_name: Name of the program (e.g., 'phenix.refine')

    Returns:
        ProgramValidator: The validator, or None if not found

    Example:
        validator = get_validator('phenix.refine')
        if validator:
            command = validator.prevalidate(command, state)
    """
    return _VALIDATORS.get(program_name)


def get_all_validators() -> dict:
    """
    Get all registered validators.

    Returns:
        dict: Mapping of program_name -> validator
    """
    return _VALIDATORS.copy()


def register_validator(validator: ProgramValidator) -> None:
    """
    Register a new validator.

    Args:
        validator: The validator instance to register

    Example:
        class MyValidator(ProgramValidator):
            @property
            def program_name(self):
                return "phenix.my_program"

        register_validator(MyValidator())
    """
    _VALIDATORS[validator.program_name] = validator
