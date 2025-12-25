"""
Document loading and chunking for RAG database building.

This module handles:
- Loading documents from folders (TXT, PDF, HTML)
- Custom HTML loader for Phenix documentation
- Smart chunking with keyword section handling

Usage:
    from libtbx.langchain.rag import load_all_docs_from_folder

    docs, files = load_all_docs_from_folder('./docs/')
"""
from __future__ import absolute_import, division, print_function

import os
import re
from typing import List

from bs4 import BeautifulSoup
from langchain_core.documents import Document
from langchain_community.document_loaders import TextLoader, PyPDFLoader, UnstructuredHTMLLoader
from langchain_community.document_loaders.base import BaseLoader
from langchain_text_splitters import RecursiveCharacterTextSplitter


# =============================================================================
# Custom Loaders
# =============================================================================

class PhenixHTMLLoader(BaseLoader):
    """
    Custom loader that parses HTML and ensures newlines are inserted
    between tags (like list items), preventing 'one long line' issues.
    """
    def __init__(self, file_path: str):
        self.file_path = file_path

    def load(self) -> List[Document]:
        with open(self.file_path, "r", encoding="utf-8", errors="ignore") as f:
            soup = BeautifulSoup(f, "html.parser")

        for script in soup(["script", "style"]):
            script.extract()

        text = soup.get_text(separator="\n")
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        cleaned_text = "\n".join(lines)

        metadata = {"source": self.file_path}
        return [Document(page_content=cleaned_text, metadata=metadata)]


# =============================================================================
# Chunking Functions
# =============================================================================

def _custom_chunker(docs: List[Document], keyword_phrase: str = "List of all available keywords") -> List[Document]:
    """
    Chunks documents with UNIVERSAL context injection.
    Prepends the source filename to EVERY chunk.

    Args:
        docs: List of documents to chunk
        keyword_phrase: Phrase that triggers special keyword section handling

    Returns:
        List of chunked documents with source context headers
    """
    final_chunks = []

    standard_splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200)
    keyword_splitter = RecursiveCharacterTextSplitter(chunk_size=2000, chunk_overlap=200)

    for doc in docs:
        content = doc.page_content
        metadata = doc.metadata

        source_path = metadata.get('source', '')
        filename = os.path.basename(source_path)
        context_header = f"** DOCUMENT: {filename} **\n"

        match = re.search(r"list\s+of\s+all\s+available\s+keywords", content, re.IGNORECASE)
        start_index = match.start() if match else -1

        current_doc_chunks = []

        if start_index != -1:
            # Part 1: Intro
            before_content = content[:start_index]
            if before_content.strip():
                before_chunks = standard_splitter.create_documents(
                    [before_content], metadatas=[metadata])
                current_doc_chunks.extend(before_chunks)

            # Part 2: Keywords
            keyword_chunk_content = content[start_index:]
            kw_chunks = keyword_splitter.create_documents(
                [keyword_chunk_content], metadatas=[metadata])
            current_doc_chunks.extend(kw_chunks)

        else:
            # Part 1: All content
            current_doc_chunks = standard_splitter.split_documents([doc])

        # --- UNIVERSAL STAMPING ---
        for chunk in current_doc_chunks:
            chunk.page_content = context_header + chunk.page_content
            final_chunks.append(chunk)

    return final_chunks


# =============================================================================
# Document Loading Functions
# =============================================================================

def load_all_docs_from_folder(folder_path: str,
                              excluded_dirs: list = None) -> tuple:
    """
    Loads all supported documents from a folder recursively.

    Args:
        folder_path: Path to folder containing documents
        excluded_dirs: List of directory names to skip

    Returns:
        tuple: (list of Documents, list of processed file paths)

    Example:
        docs, files = load_all_docs_from_folder('./phenix_docs/', excluded_dirs=['old', 'archive'])
        print(f"Loaded {len(docs)} documents from {len(files)} files")
    """
    print(f"Loading documents from {folder_path}...")

    if excluded_dirs is None:
        excluded_dirs = set()
    else:
        excluded_dirs = set(excluded_dirs)

    all_docs = []
    processed_files = []
    loaders = {
        '.txt': TextLoader,
        '.pdf': PyPDFLoader,
        '.html': PhenixHTMLLoader,
    }

    for dirpath, dirnames, filenames in os.walk(folder_path, topdown=True):
        dirnames[:] = [d for d in dirnames if d not in excluded_dirs]
        for filename in filenames:
            if filename.startswith('.'):
                continue

            file_path = os.path.join(dirpath, filename)
            file_ext = os.path.splitext(filename)[1].lower()

            if file_ext in loaders:
                loader_cls = loaders[file_ext]
                try:
                    loader = loader_cls(file_path)
                    all_docs.extend(loader.load())
                    processed_files.append(file_path)
                    print("LOADING: %s (%s)" % (file_path, len(all_docs)))
                except Exception as e:
                    print(f"Error loading file {file_path}: {e}")

    print(f"Loaded {len(all_docs)} documents.")
    return all_docs, processed_files


def load_specific_docs(file_path_list: List[str]) -> List[Document]:
    """
    Loads specific documents from a list of file paths.

    Args:
        file_path_list: List of file paths to load

    Returns:
        List of loaded Documents

    Example:
        docs = load_specific_docs(['doc1.pdf', 'doc2.html'])
    """
    all_docs = []
    loaders = {
        '.txt': TextLoader,
        '.pdf': PyPDFLoader,
        '.html': UnstructuredHTMLLoader
    }

    print(f"Loading {len(file_path_list)} specific files...")

    for file_path in file_path_list:
        file_ext = os.path.splitext(file_path)[1].lower()
        if file_ext in loaders:
            loader_cls = loaders[file_ext]
            try:
                loader = loader_cls(file_path)
                all_docs.extend(loader.load())
                print(f"  - Successfully loaded: {os.path.basename(file_path)}")
            except Exception as e:
                print(f"  - WARNING: Error loading file {file_path}: {e}")
        else:
            print(f"  - WARNING: No loader available for file type: {file_path}")

    print(f"Successfully loaded content from {len(all_docs)} documents.")
    return all_docs
