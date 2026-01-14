from __future__ import absolute_import, division, print_function
from langgraph.graph import StateGraph, END

# Import from libtbx.langchain
from libtbx.langchain.agent.graph_state import AgentState
from libtbx.langchain.agent.graph_nodes import (
    perceive, plan, build, validate, fallback, output_node
)


def route_after_validate(state):
  """
  Conditional Edge Logic after Validate node.

  Routes to:
  - "output": if stop flag set or validation passed
  - "plan": if validation failed and attempts < 3 (retry)
  - "fallback": if validation failed and attempts >= 3
  """
  # Check stop flag first (workflow complete)
  if state.get("stop"):
    return "output"

  # Check validation result
  if state.get("validation_error"):
    attempt = state.get("attempt_number", 0)
    if attempt < 3:
      return "plan"  # Retry with feedback
    else:
      return "fallback"  # Give up, use mechanical

  # Validation passed
  return "output"


def build_agent_graph():
  """
  Build and compile the LangGraph agent.

  Graph topology:

    perceive --> plan --> build --> validate --+--> output --> END
                  ^                            |
                  |  (retry if attempts < 3)   |
                  +----------------------------+
                                               |
                              (fallback if attempts >= 3)
                                               |
                                               v
                                           fallback --> output --> END

  Returns:
    Compiled LangGraph application
  """
  workflow = StateGraph(AgentState)

  # Add Nodes
  workflow.add_node("perceive", perceive)
  workflow.add_node("plan", plan)
  workflow.add_node("build", build)
  workflow.add_node("validate", validate)
  workflow.add_node("fallback", fallback)
  workflow.add_node("output", output_node)

  # Linear Edges
  workflow.set_entry_point("perceive")
  workflow.add_edge("perceive", "plan")
  workflow.add_edge("plan", "build")
  workflow.add_edge("build", "validate")

  # Conditional Edge from Validate
  workflow.add_conditional_edges(
    "validate",
    route_after_validate,
    {
      "plan": "plan",
      "fallback": "fallback",
      "output": "output"
    }
  )

  # Fallback always goes to output
  workflow.add_edge("fallback", "output")

  # Output goes to END
  workflow.add_edge("output", END)

  return workflow.compile()


# Convenience function for external use
def invoke_agent(available_files, log_text="", history=None, user_advice="", max_cycles=20):
  """
  Convenience function to invoke the agent graph.

  Args:
    available_files: List of files available on client
    log_text: Log text to analyze
    history: List of previous cycle records
    user_advice: User instructions
    max_cycles: Maximum cycles

  Returns:
    dict: Final state with command and reasoning
  """
  from graph_state import create_initial_state

  app = build_agent_graph()

  initial_state = create_initial_state(
    available_files=available_files,
    log_text=log_text,
    history=history,
    user_advice=user_advice,
    max_cycles=max_cycles
  )

  return app.invoke(initial_state)
