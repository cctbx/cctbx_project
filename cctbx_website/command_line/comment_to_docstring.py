import ast
import io
import sys
from typing import List, Tuple, Optional

def is_special_comment(line: str) -> bool:
    """Checks if a line is a special comment (shebang or coding) that should be ignored."""
    stripped_line = line.strip()
    if stripped_line.startswith('#!'):
        return True
    if 'coding:' in stripped_line:
        return True
    return False

def find_module_comment_block(lines: List[str]) -> Optional[Tuple[str, int, int]]:
    """Finds the first top-level (unindented) comment block in the file."""
    first_comment_index = -1
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith(('class ', 'def ')):
            break
        if is_special_comment(line) or not stripped:
            continue
        if stripped.startswith('#') and (len(line) - len(line.lstrip(' '))) == 0:
            if first_comment_index == -1:
                first_comment_index = i
    if first_comment_index == -1:
        return None

    end_of_block_index = first_comment_index
    for i in range(first_comment_index + 1, len(lines)):
        line = lines[i]
        stripped_line = line.strip()
        if stripped_line and not stripped_line.startswith('#'):
            break
        if stripped_line.startswith('#') and (len(line) - len(line.lstrip(' '))) != 0:
            break
        end_of_block_index = i

    content = " ".join([lines[i].strip().lstrip('#').strip() for i in range(first_comment_index, end_of_block_index + 1) if lines[i].strip().startswith('#')])
    return content, first_comment_index, end_of_block_index

def find_local_comment_block(lines: List[str], start_index: int, end_index: int, indent_level: int) -> Optional[Tuple[str, int, int]]:
    """Finds a comment block for a function/class within a specific range and indentation."""
    first_comment_line = -1
    for i in range(start_index, end_index):
        line = lines[i]
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith('#'):
            if len(line) - len(line.lstrip(' ')) == indent_level:
                first_comment_line = i
                break
            else:
                return None
        else:
            return None
    if first_comment_line == -1:
        return None

    last_comment_line = first_comment_line
    for i in range(first_comment_line + 1, end_index):
        line = lines[i]
        stripped = line.strip()
        if not stripped:
            last_comment_line = i
            continue
        if stripped.startswith('#') and len(line) - len(line.lstrip(' ')) == indent_level:
            last_comment_line = i
            continue
        break

    content_list = [lines[i].strip().lstrip('#').strip() for i in range(first_comment_line, last_comment_line + 1) if lines[i].strip().startswith('#')]
    return " ".join(content_list), first_comment_line, last_comment_line

def convert_comments_to_docstrings(source_code: str,
      remove_commented_import: bool = True) -> str:
    """Main conversion function to process source code and convert comments to docstrings."""
    source_code = source_code.replace('\u00A0', ' ')

    if remove_commented_import:
      """Replace any comment lines that start with '#' in column 1 and contain
      the word "import " with blank lines.  These are not interesting and
      otherwise show up as doc strings. """
      lines = source_code.splitlines()
      for i in range(len(lines)):
        if lines[i].startswith("#") and lines[i].find("import ") > -1:
          lines[i] = ""
      source_code = "\n".join(lines)

    try:
        tree = ast.parse(source_code)
    except SyntaxError as e:
        print(f"Error: Could not parse source. Invalid syntax on line {e.lineno}.", file=sys.stderr)
        return source_code

    lines = source_code.splitlines()
    original_had_trailing_newline = source_code.endswith('\n') or source_code.endswith('\r\n')
    replacements = []

    if not ast.get_docstring(tree):
        result = find_module_comment_block(lines)
        if result:
            replacements.append({"type": "module", "content": result[0], "start_idx": result[1], "end_idx": result[2]})

    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.ClassDef)) or ast.get_docstring(node) or not node.body:
            continue

        signature_end_lineno = -1
        for i in range(node.lineno - 1, len(lines)):
            if lines[i].strip().endswith(':'):
                signature_end_lineno = i + 1
                break
        if signature_end_lineno == -1:
            continue

        first_body_line = lines[node.body[0].lineno - 1]
        expected_indent = len(first_body_line) - len(first_body_line.lstrip(' '))

        search_start_index = signature_end_lineno
        search_end_index = node.body[0].lineno - 1

        if search_start_index > search_end_index:
            continue

        result = find_local_comment_block(lines, search_start_index, search_end_index, expected_indent)
        if result:
            replacements.append({"type": "local", "content": result[0], "start_idx": result[1], "end_idx": result[2], "indent": " " * expected_indent, "name": f"{node.__class__.__name__} '{node.name}'"})

    module_rep = next((r for r in replacements if r["type"] == "module"), None)
    if module_rep:
        start_idx, end_idx, content = module_rep['start_idx'], module_rep['end_idx'], module_rep['content']

        del lines[start_idx : end_idx + 1]

        num_lines = end_idx - start_idx + 1
        padded_content = f" {content} "
        docstring_lines = [f'"""{padded_content}"""'] if num_lines == 1 else [f'"""{padded_content.rstrip()}'] + [''] * (num_lines - 2) + ['"""']

        # --- FIX: New robust logic for finding the insertion point ---
        insert_pos = 0
        for i, line in enumerate(lines):
            # Keep advancing past special comments or blank lines at the top of the file
            if is_special_comment(line) or not line.strip():
                insert_pos = i + 1
            else:
                # Stop at the first line of actual code
                break

        lines[insert_pos:insert_pos] = docstring_lines

    local_reps = sorted([r for r in replacements if r['type'] == 'local'], key=lambda x: x['start_idx'], reverse=True)
    for rep in local_reps:
        start_idx, end_idx, content, indent, name = rep['start_idx'], rep['end_idx'], rep['content'], rep['indent'], rep['name']
        num_lines = end_idx - start_idx + 1

        if num_lines == 1:
            content_for_doc = f" {content} "
            new_lines = [f'{indent}"""{content_for_doc}"""']
        else:
            content_for_doc = f" {content}"
            new_lines = [f'{indent}"""{content_for_doc}']
            new_lines.extend([''] * (num_lines - 2))
            new_lines.append(f'{indent}"""')

        lines[start_idx : end_idx+1] = new_lines

    output = "\n".join(lines)

    if original_had_trailing_newline and len(output) > 0 and not output.endswith('\n'):
        return output + '\n'

    return output

if __name__ == '__main__':
  args = sys.argv[1:]
  if 'dry_run' in args:
     dry_run = True
     args.remove('dry_run')
  else:
     dry_run = False

  filepath = sys.argv[1]
  with open(filepath, 'r', encoding='utf-8') as f:
    original_code = f.read()
  modified_code = convert_comments_to_docstrings(original_code)
  if not dry_run:
    print(modified_code)

