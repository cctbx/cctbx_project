"""
Core type definitions for the Phenix Crystallography Agent.

These dataclasses define the data structures used throughout the system.
Using dataclasses provides:
- Automatic __init__, __repr__, __eq__
- Type hints for IDE support
- Clean, readable code
"""
from __future__ import absolute_import, division, print_function
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class AgentPlan:
    """
    A plan for what action(s) to take next.

    Attributes:
        commands: List of Phenix commands to execute
        reasoning: Explanation of why this plan was chosen
        program: The primary program being run (e.g., 'phenix.refine')
        strategy: Strategy details for command construction
        stop_conditions: Conditions that should stop execution
        estimated_cycles: Expected number of cycles to complete

    Example:
        plan = AgentPlan(
            commands=["phenix.refine model.pdb data.mtz"],
            reasoning="Model needs refinement after MR",
            program="phenix.refine",
            strategy="initial refinement with default parameters"
        )
    """
    commands: list = field(default_factory=list)
    reasoning: str = ""
    program: str = ""
    strategy: str = ""
    stop_conditions: list = field(default_factory=list)
    estimated_cycles: int = 1


@dataclass
class ValidationResult:
    """
    Result of command validation.

    Attributes:
        is_valid: Whether the command passed validation
        error_message: Description of validation error (if any)
        original_command: The command before any fixes
        fixed_command: The command after fixes (if any were applied)
        fix_attempts: Number of attempts made to fix the command

    Example:
        result = ValidationResult(
            is_valid=True,
            original_command="phenix.refine model.pdb data.mtz",
            fixed_command="phenix.refine model.pdb data.mtz"
        )
    """
    is_valid: bool = False
    error_message: str = ""
    original_command: str = ""
    fixed_command: str = ""
    fix_attempts: int = 0


@dataclass
class CommandResult:
    """
    Result of executing a Phenix command.

    Attributes:
        command: The command that was executed
        program: The program name (e.g., 'phenix.refine')
        success: Whether the command completed without errors
        error: Error message if command failed
        log_file: Path to the log file
        output_files: List of files produced by the command
        summary: Summary of what the command accomplished
        analysis: Detailed analysis of the results

    Example:
        result = CommandResult(
            command="phenix.refine model.pdb data.mtz",
            program="phenix.refine",
            success=True,
            log_file="/path/to/refine.log",
            output_files=["model_refine_001.pdb", "model_refine_001.mtz"]
        )
    """
    command: str = ""
    program: str = ""
    success: bool = False
    error: str = ""
    log_file: str = ""
    output_files: list = field(default_factory=list)
    summary: str = ""
    analysis: str = ""


@dataclass
class ProjectState:
    """
    Current state of a crystallography project.

    Tracks all known information about the project including
    data files, models, and crystallographic parameters.

    Attributes:
        crystallography_or_cryoem: Type of project ('crystallography' or 'cryoem')
        sequence_file: Path to sequence file
        crystallography: Dict of crystallographic parameters
        cryoem: Dict of cryo-EM parameters
        model_list: List of model files with metadata

    Example:
        state = ProjectState(
            crystallography_or_cryoem='crystallography',
            sequence_file='/path/to/seq.fa',
            crystallography={
                'space_group': 'P 21 21 21',
                'resolution': 2.0,
                'original_data': '/path/to/data.mtz'
            },
            model_list=[
                {'file': 'model.pdb', 'r_work': 0.25, 'r_free': 0.28}
            ]
        )
    """
    crystallography_or_cryoem: str = "crystallography"
    sequence_file: str = ""
    crystallography: dict = field(default_factory=dict)
    cryoem: dict = field(default_factory=dict)
    model_list: list = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict) -> 'ProjectState':
        """Create ProjectState from a dictionary."""
        return cls(
            crystallography_or_cryoem=data.get('crystallography_or_cryoem', 'crystallography'),
            sequence_file=data.get('sequence_file', ''),
            crystallography=data.get('crystallography', {}),
            cryoem=data.get('cryoem', {}),
            model_list=data.get('model_list', [])
        )

    def to_dict(self) -> dict:
        """Convert ProjectState to a dictionary."""
        return {
            'crystallography_or_cryoem': self.crystallography_or_cryoem,
            'sequence_file': self.sequence_file,
            'crystallography': self.crystallography,
            'cryoem': self.cryoem,
            'model_list': self.model_list
        }


@dataclass
class NextMove:
    """
    Complete result from the agent's planning process.

    This is the main output of generate_next_move(), containing
    the command to run and all supporting information.

    Attributes:
        command: The Phenix command to execute
        program: The program name
        explanation: Full explanation including plan and reasoning
        strategy: Strategy details
        process_log: Log of the agent's thought process
        error: Error message if planning failed

    Example:
        move = NextMove(
            command="phenix.refine model.pdb data.mtz",
            program="phenix.refine",
            explanation="PLAN: Refine model\\nREASONING: Improve R-factors",
            strategy="default refinement"
        )
    """
    command: str = ""
    program: str = ""
    explanation: str = ""
    strategy: str = ""
    process_log: str = ""
    error: Optional[str] = None

    @property
    def success(self) -> bool:
        """Returns True if planning succeeded (no error)."""
        return self.error is None and self.command != ""

