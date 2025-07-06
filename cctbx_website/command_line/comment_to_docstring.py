from __future__ import absolute_import, division, print_function
import ast
import sys
from typing import List, Tuple, Optional

def find_description_in_text(source_code):
    import ast
    import sys

    try:
        tree = ast.parse(source_code)
    except SyntaxError as e:
        return ""

    # Visitor class to traverse the Abstract Syntax Tree
    class DescriptionVisitor(ast.NodeVisitor):
        def __init__(self):
            self.descriptions_found = []

        def visit_Assign(self, node):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == 'description':
                    # For Python 3.8+, string literals are ast.Constant
                    if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                        self.descriptions_found.append(node.value.value)
                    # For Python 3.0-3.7, they are ast.Str
                    elif isinstance(node.value, ast.Str):
                        self.descriptions_found.append(node.value.s)
            # Continue traversing to find all instances in the file
            self.generic_visit(node)

    # Create an instance of the visitor and run it on the AST
    visitor = DescriptionVisitor()
    visitor.visit(tree)

    # Join the collected descriptions into a single string
    return "\n".join(visitor.descriptions_found)


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

def remove_commented_imports(source_code):
      """Replace any comment lines that start with '#' in column 1 and contain
      the word "import " or other indications of non-useful information
      with blank lines.  These are not interesting and
      otherwise show up as doc strings. """
      non_useful_text = ["import ", "LIBTBX_","BOOST_ADAPTBX_","=========",
         "-----------","__________"]
      lines = source_code.splitlines()
      for i in range(len(lines)):
        if lines[i].startswith("#"):
          for x in non_useful_text:
            if lines[i].find(x) > -1:
              lines[i] = ""
              break
      source_code = "\n".join(lines)
      return source_code

def add_description_if_available(source_code):
  """ If a description variable is set with text in the source code,
  and no docstring is present for the file as a whole,
  put in a docstring equal to that description at the top"""

  if not source_code:
    return source_code

  ss = source_code.strip()
  if ss.startswith('"""') or ss.startswith("'''"):
    return source_code
  desc = find_description_in_text(source_code)
  if desc:
    source_code = '"""%s"""\n%s' %(desc, source_code)
  return source_code


def convert_comments_to_docstrings(source_code: str,
      remove_commented_import: bool = True,
      raise_on_errors: bool = False) -> str:
    """Main conversion function to process source code and
        convert comments to docstrings."""

    # Clean up to start
    source_code = source_code.replace('\u00A0', ' ')

    if remove_commented_import:
      source_code = remove_commented_imports(source_code)
      source_code = add_description_if_available(source_code)


    try:
        tree = ast.parse(source_code)
    except SyntaxError as e:
        print(f"Error: Could not parse source. Invalid syntax on line {e.lineno}.", file=sys.stderr)
        if raise_on_errors:
          raise AssertionError("Failed to parse source")
        else: # usual
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
