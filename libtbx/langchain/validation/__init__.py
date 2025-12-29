"""
Validation module - command validation and fixing.

This module contains:
- base.py: Abstract base class for validators
- core_validator.py: Generic validation functions
- mtz_utils.py: MTZ file utilities
- phenix_refine.py: phenix.refine-specific validator
- registry.py: Validator lookup and registration
"""
from __future__ import absolute_import, division, print_function

from libtbx.langchain.validation.base import ProgramValidator

from libtbx.langchain.validation.core_validator import (
    validate_phenix_command,
    fix_command_syntax,
    validate_is_phenix_command,
)

from libtbx.langchain.validation.mtz_utils import (
    mtz_has_rfree_flags,
    add_rfree_generation_if_needed,
)

from libtbx.langchain.validation.phenix_refine import RefineValidator

from libtbx.langchain.validation.registry import (
    get_validator,
    get_all_validators,
    register_validator,
)

__all__ = [
    # Base
    'ProgramValidator',
    # Core validation
    'validate_phenix_command',
    'fix_command_syntax',
    'validate_is_phenix_command',
    # MTZ utilities
    'mtz_has_rfree_flags',
    'add_rfree_generation_if_needed',
    # Validators
    'RefineValidator',
    # Registry
    'get_validator',
    'get_all_validators',
    'register_validator',
]
