"""
Base class for planning strategies.

A strategy decides what Phenix program(s) to run next based on:
- Run history (what has been done)
- Project state (what we know about the structure)
- Available tools (what programs can be run)

To create a new strategy:
1. Subclass PlanningStrategy
2. Implement the abstract methods
3. Register in strategies/__init__.py

Example:
    class MyStrategy(PlanningStrategy):
        @property
        def name(self) -> str:
            return "my_strategy"

        @property
        def description(self) -> str:
            return "My custom planning approach"

        async def plan_next_move(self, run_history, project_state, context):
            # Your planning logic here
            return AgentPlan(
                commands=["phenix.refine model.pdb data.mtz"],
                reasoning="Need to refine the model",
                program="phenix.refine"
            )
"""
from __future__ import absolute_import, division, print_function

from abc import ABC, abstractmethod
from typing import Any

from libtbx.langchain.core.types import AgentPlan


class PlanningStrategy(ABC):
    """
    Abstract base class for all planning strategies.

    Subclasses must implement:
    - name: Unique identifier for this strategy
    - plan_next_move(): Generate the next action plan

    Optional overrides:
    - description: Human-readable description
    - initialize(): Setup called once before first use
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """
        Unique name for this strategy.

        Used for:
        - Strategy selection: CrystallographyAgent(strategy='sequential')
        - Logging and debugging
        - Registry lookup

        Returns:
            str: Short identifier (e.g., 'sequential', 'batch', 'workflow')
        """
        pass

    @property
    def description(self) -> str:
        """
        Human-readable description of what this strategy does.

        Returns:
            str: Description for documentation/help
        """
        return f"Planning strategy: {self.name}"

    def initialize(self, llm: Any, embeddings: Any, db_dir: str) -> None:
        """
        Initialize the strategy with required resources.

        Called once before first use. Override to set up any
        resources needed for planning (e.g., RAG retriever).

        Args:
            llm: Language model for reasoning
            embeddings: Embedding model for RAG
            db_dir: Path to documentation database
        """
        self.llm = llm
        self.embeddings = embeddings
        self.db_dir = db_dir

    @abstractmethod
    async def plan_next_move(
        self,
        run_history: list,
        project_state: dict,
        context: dict = None
    ) -> AgentPlan:
        """
        Generate a plan for what to do next.

        This is the main method that each strategy must implement.

        Args:
            run_history: List of previous runs with their results
            project_state: Current state of the project (files, parameters)
            context: Additional context (advice, original_files, etc.)

        Returns:
            AgentPlan: The planned action(s) to take

        Example implementation:
            async def plan_next_move(self, run_history, project_state, context):
                # Analyze what's been done
                if not run_history:
                    return AgentPlan(
                        commands=["phenix.xtriage data.mtz"],
                        reasoning="First step: analyze data quality",
                        program="phenix.xtriage"
                    )

                # Plan based on history
                last_run = run_history[-1]
                if last_run.get('program') == 'phenix.xtriage':
                    return AgentPlan(
                        commands=["phenix.phaser ..."],
                        reasoning="Data looks good, try MR",
                        program="phenix.phaser"
                    )
        """
        pass

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}(name='{self.name}')>"

