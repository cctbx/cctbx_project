"""
word_index_generator (generate an index from html directories)

Generate a multi-page alphabetical HTML word index from a directory tree of HTML files.
Only visible, non-trivial words (excluding stopwords and filtered tokens) are indexed.
Each section page contains navigable links and structured grouping by word prefix.
Formatted to match pdoc3-style HTML structure.
"""
from __future__ import division

import os
import re
import html

OK = True
try:
  from bs4 import BeautifulSoup, Comment
  from collections import defaultdict
  from nltk.corpus import stopwords
  import nltk
except Exception as e:
  OK = False

def download_stopwords():
    """Ensure stopwords are available."""
    nltk.download('stopwords')

def get_visible_text(html_content):
    """
    Extracts visible text from HTML content, excluding scripts, styles, comments, and hidden elements.
    """
    soup = BeautifulSoup(html_content, "html.parser")
    for tag in soup(['pre', 'script', 'style', 'meta', 'noscript']):
        tag.decompose()
    for comment in soup.find_all(string=lambda text: isinstance(text, Comment)):
        comment.extract()
    for tag in soup.find_all(style=True):
        style = tag['style'].replace(" ", "").lower()
        if "display:none" in style or "visibility:hidden" in style:
            tag.decompose()
    body = soup.body or soup
    return body.get_text(separator=' ', strip=True)

def tokenize(text):
    """Tokenize input text into lowercase words."""
    return re.findall(r'\b\w+\b', text.lower())

def get_letter_section(word):
    """Returns the first letter of the word or 'other' if not alphabetic."""
    first = word[0].lower()
    return first if first.isalpha() else "other"

def build_word_index(html_dir, stop_words, exclude_pattern, index_dir):
    """
    Walk the HTML directory and build a word index {word -> set of file paths (up to 5)}
    """
    word_index = defaultdict(set)
    for root, _, files in os.walk(html_dir):
        for file in files:
            if file.endswith(".html"):
                filepath = os.path.abspath(os.path.join(root, file))
                rel_filepath = os.path.join("..", os.path.relpath(filepath, os.path.abspath(html_dir)))
                if filepath.find(index_dir) > -1:
                    continue
                try:
                    with open(filepath, "r", encoding="utf-8") as f:
                        html_content = f.read()
                        text = get_visible_text(html_content)

                        if not text.strip():
                            print(f"⚠️ No visible text in: {filepath}")
                            continue
                        print(f"✅ Visible text found in: {filepath}")

                        words = tokenize(text)
                        for word in words:
                            if (
                                word not in stop_words and
                                len(word) > 2 and
                                not word[0].isdigit() and
                                not word.startswith("_") and
                                not exclude_pattern.match(word) and
                                word.isascii()
                            ):
                                if len(word_index[word]) < 5:
                                    word_index[word].add(rel_filepath)
                except Exception as e:
                    print(f"⚠️ Skipped {filepath}: {e}")
    return word_index

def generate_html_pages(word_index, index_title, index_dir):
    """
    Generate alphabetically sectioned HTML index files with navigation and pdoc3-style layout.
    """
    sorted_words = sorted(word_index.items())
    sectioned_words = defaultdict(list)
    for word, paths in sorted_words:
        section = get_letter_section(word)
        sectioned_words[section].append((word, paths))

    sections_sorted = sorted(sectioned_words.keys())
    overall_file = None
    for section, word_list in sectioned_words.items():
        section_file = os.path.join(index_dir, f"word_index_{section}.html")
        current_index = sections_sorted.index(section)
        prev_section = sections_sorted[current_index - 1] if current_index > 0 else None
        next_section = sections_sorted[current_index + 1] if current_index + 1 < len(sections_sorted) else None

        lines = [
            "<!doctype html>",
            "<html lang=\"en\">",
            "<head>",
            "<meta charset=\"utf-8\">",
            f"<title>{index_title} Index: {section.upper()}</title>",
            "<meta name='viewport' content='width=device-width, initial-scale=1'>",
            "<link rel=\"stylesheet\" href=\"https://cdnjs.cloudflare.com/ajax/libs/10up-sanitize.css/13.0.0/sanitize.min.css\">",
            "<link rel=\"stylesheet\" href=\"https://cdnjs.cloudflare.com/ajax/libs/10up-sanitize.css/13.0.0/typography.min.css\">",
            "<style>",
            "body { font-family: system-ui, sans-serif; background-color: #fff; color: #222; line-height: 1.6; margin: 0; }",
            "main { display: flex; flex-direction: row; flex-wrap: wrap; }",
            "#sidebar { width: 15%; min-width: 160px; padding: 1.5em; background-color: #f9f9f9; border-right: 1px solid #ddd; position: sticky; top: 0; height: 100vh; overflow-y: auto; }",
            "article#content { flex: 1; min-width: 300px; padding: 3em 4em; }",
            "#jump-list { column-count: 2; -webkit-column-count: 2; -moz-column-count: 2; padding-left: 0; list-style: none; }",
            "article ul { column-count: 3; -webkit-column-count: 3; -moz-column-count: 3; padding-left: 0; list-style: none; }",
            "li { margin-bottom: .5em; }",
            "h1, h2, h3 { font-weight: 300; color: #111; }",
            "code { font-family: 'DejaVu Sans Mono', monospace; background-color: #f3f3f3; padding: 1px 4px; border-radius: 3px; }",
            ".nav-links { margin-bottom: 1em; font-size: 0.9em; }",
            "footer { font-size: 0.75em; padding: 1em; text-align: right; color: #999; }",
            "@media (max-width: 768px) { #sidebar { position: relative; height: auto; } article#content { padding: 1em; } }",

            "a { color: #336699; /* softer blue-gray tint */ text-decoration: none; }",
            "a:hover { color: #224466; /* darker shade on hover */ }",

            "</style>",
            "</head>",
            "<body>",
            "<main>",
            "<nav id=\"sidebar\">",
            "<ul id=\"index\">",
            f"<li><a href=\"../index.html\">{index_title}</a></li>",
            "</ul></li>",
            "<li><h3>Jump to</h3>",
            "<ul id=\"jump-list\">",
        ]

        for sec in sections_sorted:
            lines.append(f"<li><a href='word_index_{sec}.html'>{sec.upper()}</a></li>")

        lines += [
            "</ul></li>",
            "</ul>",
            "</nav>",
            "<article id=\"content\">",
            f"<h1>{index_title} Index: {section.upper()}</h1>",
            "<div class=\"nav-links\">",
        ]

        if prev_section:
            lines.append(f"<a href='word_index_{prev_section}.html'>&larr; Previous ({prev_section.upper()})</a>")
        if next_section:
            lines.append(f" | <a href='word_index_{next_section}.html'>Next ({next_section.upper()}) &rarr;</a>")

        lines.append("</div>")

        lines.append("<div class='nav-links'><strong>Subsections:</strong> ")
        subsection_prefixes = sorted(set(word[:2].lower() for word, _ in word_list))
        for prefix in subsection_prefixes:
            lines.append(f"<a href='#{prefix}'>{prefix}</a> ")
        lines.append("</div>")

        current_prefix = ""
        first_prefix = True

        for word, paths in word_list:
            prefix = word[:2].lower()
            if prefix != current_prefix:
                if not first_prefix:
                    lines.append("</ul>")
                lines.append(f"<h3 id='{prefix}'>{prefix}</h3><ul>")
                current_prefix = prefix
                first_prefix = False

            word_html = f"<strong><code>{html.escape(word)}</code></strong>: " + ", ".join(
                f"<a href='{html.escape(path)}' target='_blank'>{os.path.basename(path)}</a>"
                for path in sorted(paths)
            )
            lines.append(f"<li>{word_html}</li>")

        lines.append("</ul>")
        lines.append("</article></main>")
        lines.append("<footer id=\"footer\"><p>Generated by word_index_generator</p></footer>")
        lines.append("</body></html>")

        with open(section_file, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))

        if section == "a":
            overall_file = os.path.join(index_dir, "index.html")
            with open(overall_file, "w", encoding="utf-8") as f:
                f.write("\n".join(lines))

        print(f"✅ Created section: {section_file}")

    print(f"✅ Main index written to '{overall_file}'")

def run(args):
    """
    Entry point: index HTML files in a directory and generate navigable A–Z word index.
    Expects path to directory containing html files, path for index files
    (normally same as path to html files + /index_files/ and
    optional third arg with title.
    """
    try:
        html_dir = args[0]
        assert os.path.isdir(html_dir)
        index_dir = args[1]
        assert os.path.relpath(index_dir, html_dir).find(os.path.sep) < 0
        if not os.path.isdir(index_dir):
            os.mkdir(index_dir)
        print(f"HTML to be read from '{html_dir}'")
        print(f"Indexing HTML to be written to '{index_dir}'")
    except Exception as e:
        print("Please run with path to directory containing html files as 1st arg and path to directory inside that for index files (index_files usually)")
        return

    index_title = args[2] if len(args) > 2 else "Word"
    print(f"Title to use: '{index_title}'")

    download_stopwords()
    stop_words = set(stopwords.words('english'))
    exclude_pattern = re.compile(r'^[a-zA-Z][0-9_]')

    word_index = build_word_index(html_dir, stop_words, exclude_pattern, index_dir)
    print(f"\n📚 Indexed {len(word_index)} unique words.")
    generate_html_pages(word_index, index_title, index_dir)

if __name__ == "__main__":
    if OK:
        import sys
        run(sys.argv[1:])
    else:
        print("Cannot run word_index_generator without nltk beautifulsoup4")
