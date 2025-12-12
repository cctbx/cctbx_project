"""
Text processing utilities.

This module contains helper functions for:
- Extracting text blocks from log files
- Processing log summaries
"""
from __future__ import absolute_import, division, print_function


def find_text_block(log_text: str, target_text: str, end_text: str = "**") -> str:
    """
    Extracts a block of text from a log starting at target_text and ending at end_text.

    Args:
        log_text: The full log text
        target_text: The text that marks the start of the block
        end_text: The text that marks the end of the block

    Returns:
        str: The extracted text block (lines joined with spaces)

    Example:
        log = '''
        **3. Program Name**
        phenix.refine
        **4. Next Section**
        '''
        program = find_text_block(log, "**3", end_text="**")
        # Returns: "**3. Program Name** phenix.refine"
    """
    text_lines = []
    started = False
    if log_text is None:
        log_text = ""
    for line in log_text.splitlines():
        if (not started) and line.strip().replace(" ", "").startswith(target_text):
            started = True
        elif (started) and line.strip().replace(" ", "").startswith(end_text):
            break
        if started:
            text_lines.append(line)
    return " ".join(text_lines)


def get_processed_log_dict(log_text: str,
                           summary_text_block_start: str = "**8",
                           program_text_block_start: str = "**3") -> dict:
    """
    Extracts key information from a formatted log summary.

    Args:
        log_text: The log summary text
        summary_text_block_start: Marker for summary section
        program_text_block_start: Marker for program name section

    Returns:
        dict: Contains 'phenix_program' and 'summary' keys

    Example:
        info = get_processed_log_dict(summary_text)
        print(f"Program: {info['phenix_program']}")
    """
    summary = ""
    phenix_program = find_text_block(
        log_text, program_text_block_start, end_text="**")
    if not phenix_program:
        phenix_program = find_text_block(
            log_text, program_text_block_start.replace("*", "#"), end_text="##")

    dd = {
        'phenix_program': phenix_program,
        'summary': summary
    }
    return dd

