"""
Agent Session Analysis.

This module provides functions to analyze and summarize completed AI agent sessions.
It generates structured summaries and optionally includes LLM-based assessments.

Usage:
    from libtbx.langchain.analysis.agent_session_analyzer import (
        analyze_agent_session,
        generate_session_summary,
    )

    # Generate summary without LLM
    summary = generate_session_summary(session, include_llm=False)

    # Generate summary with LLM assessment
    summary = analyze_agent_session(session, llm, timeout=60)
"""
from __future__ import absolute_import, division, print_function

import asyncio
from libtbx import group_args


def generate_session_summary(session, include_llm_placeholder=False):
    """
    Generate a structured Markdown summary of an agent session.

    This function extracts data from the session and formats it as
    a clean Markdown report. No LLM is required.

    Args:
        session: AgentSession object with completed session data
        include_llm_placeholder: If True, include placeholder for LLM assessment

    Returns:
        group_args with:
            - markdown: The formatted Markdown summary
            - data: Structured data dict
            - error: Error message if any
    """
    try:
        result = session.generate_agent_session_summary(
            include_llm_assessment=include_llm_placeholder
        )

        return group_args(
            group_args_type='session_summary',
            markdown=result["markdown"],
            data=result["data"],
            error=None,
        )

    except Exception as e:
        return group_args(
            group_args_type='error',
            markdown=None,
            data=None,
            error=str(e),
        )


async def analyze_agent_session_async(session, llm, timeout=60):
    """
    Analyze an agent session with LLM assessment (async version).

    Args:
        session: AgentSession object with completed session data
        llm: Language model for assessment
        timeout: Timeout in seconds for LLM call

    Returns:
        group_args with:
            - markdown: Complete Markdown summary with assessment
            - data: Structured data dict
            - assessment: LLM-generated assessment text
            - error: Error message if any
    """
    try:
        # First generate the structured summary
        summary_result = session.generate_agent_session_summary(
            include_llm_assessment=True
        )

        markdown = summary_result["markdown"]
        data = summary_result["data"]

        # Get the concise summary for LLM
        llm_input = session.get_summary_for_llm_assessment()

        # Get the assessment prompt
        from libtbx.langchain.knowledge.prompts_hybrid import (
            get_agent_session_assessment_prompt
        )
        prompt_template = get_agent_session_assessment_prompt()
        prompt = prompt_template.format(session_summary=llm_input)

        # Call LLM for assessment
        try:
            response = await asyncio.wait_for(
                llm.ainvoke(prompt),
                timeout=timeout
            )

            if hasattr(response, 'content'):
                assessment = response.content
            else:
                assessment = str(response)

            # Replace placeholder in markdown with actual assessment
            markdown = markdown.replace(
                "_[LLM assessment will be inserted here]_",
                assessment
            )

        except asyncio.TimeoutError:
            assessment = "_Assessment timed out_"
            markdown = markdown.replace(
                "_[LLM assessment will be inserted here]_",
                assessment
            )

        except Exception as e:
            assessment = f"_Could not generate assessment: {e}_"
            markdown = markdown.replace(
                "_[LLM assessment will be inserted here]_",
                assessment
            )

        return group_args(
            group_args_type='session_analysis',
            markdown=markdown,
            data=data,
            assessment=assessment,
            error=None,
        )

    except Exception as e:
        return group_args(
            group_args_type='error',
            markdown=None,
            data=None,
            assessment=None,
            error=str(e),
        )


def analyze_agent_session(session, llm, timeout=60):
    """
    Analyze an agent session with LLM assessment (sync wrapper).

    Args:
        session: AgentSession object with completed session data
        llm: Language model for assessment
        timeout: Timeout in seconds for LLM call

    Returns:
        group_args with:
            - markdown: Complete Markdown summary with assessment
            - data: Structured data dict
            - assessment: LLM-generated assessment text
            - error: Error message if any
    """
    return asyncio.run(analyze_agent_session_async(session, llm, timeout))


def format_assessment_for_display(result):
    """
    Format the analysis result for display to the user.

    Args:
        result: group_args from analyze_agent_session

    Returns:
        str: Formatted text for display
    """
    if result.error:
        return f"Error generating session summary: {result.error}"

    lines = []
    lines.append("=" * 60)
    lines.append("AI AGENT SESSION SUMMARY")
    lines.append("=" * 60)
    lines.append("")

    if result.markdown:
        lines.append(result.markdown)

    return "\n".join(lines)
