"""
Base class for command builders.

Command builders construct Phenix commands from high-level plans.
Different builders use different approaches:
- LLMBuilder: Uses LLM to generate commands (current approach)
- TemplateBuilder: Uses predefined templates (future)

To create a new builder:
1. Subclass CommandBuilder
2. Implement the abstract methods
3. Register in commands/__init__.py

Example:
    class MyBuilder(CommandBuilder):
        @property
        def name(self) -> str:
            return "my_builder"

        async def build_command(self, program, strategy, files, state):
            # Your command building logic
            return f"{program} {files['model']} {files['data']}"
"""
from __future__ import absolute_import, division, print_function

from abc import ABC, abstractmethod
from typing import Any
assert Any is not None


class CommandBuilder(ABC):
    """
    Abstract base class for command builders.

    Subclasses must implement:
    - name: Unique identifier for this builder
    - build_command(): Construct command from plan details

    Optional overrides:
    - initialize(): Setup called once before first use
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """
        Unique name for this command builder.

        Returns:
            str: Short identifier (e.g., 'llm', 'template')
        """
        pass

    @property
    def description(self) -> str:
        """
        Human-readable description of this builder's approach.

        Returns:
            str: Description for documentation
        """
        return f"Command builder: {self.name}"

    def initialize(self, llm: Any = None, **kwargs) -> None:
        """
        Initialize the builder with required resources.

        Args:
            llm: Language model (if needed)
            **kwargs: Additional configuration
        """
        self.llm = llm

    @abstractmethod
    async def build_command(
        self,
        program: str,
        strategy: str,
        input_files: dict,
        project_state: dict,
        keywords: str = ""
    ) -> str:
        """
        Build a Phenix command from the given specifications.

        Args:
            program: Program name (e.g., 'phenix.refine')
            strategy: Strategy description (e.g., 'initial refinement')
            input_files: Dict of input files {'model': 'x.pdb', 'data': 'y.mtz'}
            project_state: Current project state
            keywords: Valid keywords for this program

        Returns:
            str: Complete command ready to execute

        Example:
            command = await builder.build_command(
                program='phenix.refine',
                strategy='refine with 5 cycles',
                input_files={'model': 'model.pdb', 'data': 'data.mtz'},
                project_state={'crystallography': {'resolution': 2.0}}
            )
            # Returns: "phenix.refine model.pdb data.mtz main.number_of_macro_cycles=5"
        """
        pass

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}(name='{self.name}')>"
