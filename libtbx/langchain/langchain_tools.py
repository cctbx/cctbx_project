"""
A library of modular functions for building and querying LangChain RAG pipelines.
Required packages: google-generativeai langchain-cohere cohere langchain langchain-google-genai langchain-chroma pypdf unstructured beautifulsoup4 tqdm markdown-it-py langchain_community linkify-it-py

This module contains methods for setting up a RAG database from files
in a directory, querying the database, summarizing a log file, and
analyzing a log file summary in the context of the database.

Prompts are supplied for queries, summarizing, and analyzing that are
specific and optimized for analyzing using the Phenix documentation


Uses asynchronous approach for sending files to Google AI.

Splits a large log file into several smaller chunks.
It then needs to ask the Google AI model to summarize each chunk.
The asyncio.gather function sends all these summarization requests to Google's API concurrently (at the same time).
The script then waits for all the responses to come back. The total wait time is roughly the time it takes for the longest single API call to complete.

"""

import nest_asyncio
nest_asyncio.apply()
import chromadb
import cohere
from cohere.core.api_error import ApiError as CohereApiError
from langchain_google_genai._common import GoogleGenerativeAIError
from langchain_openai import ChatOpenAI, OpenAIEmbeddings
import openai # For handling OpenAI-specific exceptions

import os
import time

# --- GLOBAL VARIABLE to track the last query time ---
_last_query_time = 0

os.environ['GRPC_ENABLE_FORK_SUPPORT'] = "false"   # suppress a warning message

import asyncio
from concurrent.futures import TimeoutError
from libtbx import group_args
from google.api_core import exceptions as google_exceptions
from langchain_cohere import CohereRerank
from langchain.retrievers import ContextualCompressionRetriever
from langchain_core.runnables import RunnablePassthrough
from langchain_google_genai import ChatGoogleGenerativeAI, GoogleGenerativeAIEmbeddings
from langchain_google_genai import GoogleGenerativeAIEmbeddings

from langchain_community.document_loaders import TextLoader, PyPDFLoader, UnstructuredHTMLLoader
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.chains.combine_documents import create_stuff_documents_chain
from langchain_core.prompts import PromptTemplate
from langchain_core.documents import Document
from langchain_chroma import Chroma
from typing import Iterable
from typing import List
from libtbx.langchain.run_analyze_log import save_as_html

# --- Configuration ---

# --- Prompt Management Functions ---

def get_log_map_prompt() -> PromptTemplate:
    """Returns the prompt for summarizing a single chunk of a log file."""

    template = """
    You are a highly skilled data extraction bot. Your task is to scan a chunk of a Phenix log file and pull out only the most critical information.

    **Instructions:**
    1.  **Identify Key Steps:** Look for lines that indicate the start or end of a major computational step (e.g., "Starting AlphaFold prediction", "Docking model", "Rebuilding model").
    2.  **Extract File Names:** List any input or output file names mentioned (`.pdb`, `.cif`, `.seq`, `.fasta`, `.mtz`, `.ccp4`). Ignore `.dat` files. Be sure to capture all output file names that contain the text `overall_best` as these are generally the current best result files.
    3.  **Capture Metrics:** Record any specific numbers or metrics reported, especially map-model correlation coefficients (CC), resolution estimates, scores, and R values.
    4.  **Identify X-ray vs CryoEM data** Notice whether the data are X-ray
       crystallography data (typically mtz files) or cryo-EM data (typically
       mrc or ccp4 files).
    5.  **Be Concise:** Create a brief, bulleted list of these key facts. Do not add conversational text or explanations.

    Log Chunk:
    "{text}"

    Concise, structured summary of this chunk:
    """
    return PromptTemplate.from_template(template)

def get_log_combine_prompt() -> PromptTemplate:

    """Returns the prompt for synthesizing log chunk summaries into
        a final report."""

    template = (
        "You are an expert research scientist. "
          "Synthesize the following summaries "
        "from a long log file into a single, structured report. Your report "
        "must contain all the information requested below.\n\n"
        "Summaries:\n{context}\n\n"
        "Final Report Requirements:\n"
        "1. Table of input files and their contents. Specify whether the data are X-ray or cryo-EM.\n"
        "2. Table of specified input parameters.\n"
        "3. The name of the Phenix program that was run.\n"
        "4. Detailed description of each step carried out.\n"
        "5. Detailed table of all metrics obtained.\n"
        "6. Full reproduction of any summary table found at the end of the run."
        " Do not include any information from files with the suffix of .dat\n"
        "7. List of all output files and their contents. "
         "Ignore files with the suffix of .dat.  Try to include files "
         "that contain the text `overall_best`, as these are usually "
         "the current best result files. \n"
        "8. Identification and analysis of the key output files and their "
         "evaluation metrics. Ignore files with the suffix of .dat"
    )
    return PromptTemplate(template=template, input_variables=["context"])

def get_log_analysis_prompt() -> PromptTemplate:
    """Returns the prompt for analyzing a log summary against
        the documentation."""
    template = (
        "You are expert in crystallography and cryo-EM."
        "You are a Phenix power-user. Your task is to analyze "
        "a program summary in the context of the provided documentation "
        "and research papers.\n\n"
        "Documentation Context:\n{context}\n\n"
        "Log File Summary:\n{log_summary}\n\n"
        "Based on all the information above, please perform the following "
        " analysis. Consider the events of the log summary in the broader "
        "context of a "
        "typical Phenix structure determination workflow as described in "
        "the documentation and papers, but do not describe this context. "
        "The analysis should include the following:\n\n"
        "1.For the run described in the log summary, list the inputs and "
        "briefly describe what was done."
        "Report whether the data are from crystallography (X-ray or neutron)"
        " or from cryo-EM"
        "2. List the key output files from this run, along with the values "
        "of any available metrics describing their utilities. "
        "If no metrics are available, do not provide any. "
        "3. Evaluate whether the run described in the summary was useful. "
          "List reported metrics and expected values of these metrics and "
          "consider the goals of the program.\n"
        "4. Considering whether the data are from crystallography or"
        "cryo-EM and considering the normal sequence of Phenix tool use"
        "for that type of data, suggest "
        "three concrete next steps in structure determination using Phenix "
        "or graphical tools such as Coot or Isolde that I should take, "
        "justifying each suggestion with information from the provided "
        "documentation and papers. "
        "Do not include any information on CryoFit, ShelxD, Parrot, "
        "or CCP4 tools "
        "unless there is a specific question about them. "
        "If appropriate, include validation as one of your next steps. "
        "Name the tools that are to be used, along with their inputs "
        " and outputs and what they do. "
        "Do not suggest depositing the model. "
        "Do not suggest analyzing the biological relevance. "
    )
    return PromptTemplate(template=template, input_variables=["context", "log_summary"])


def get_docs_query_prompt() -> PromptTemplate:
    """Returns the general-purpose prompt for querying the documentation RAG."""

    template = """You are an expert assistant for the Phenix software suite.
Answer the user's question based only on the following context from documentation and papers.
Provide detailed, helpful answers.
Do not discuss limitations of the sources.
Do not make up any answers.
You can say that you do not have enough information to reply.
Do not include any information on CryoFit, ShelxD, or CCP4 tools unless there
is a specific question about them.
Consider this question in the context of the process of structure determinaion in Phenix.
Focus on using Phenix tools, but include the use of Coot or Isolde if appropriate.
Name the tools that are to be used, along with their inputs and outputs and what they do.

Context:
{context}

Question:
{input}

Answer:"""
    return PromptTemplate(
        template=template, input_variables=["context", "input"]
    )


# --- Core Functionality ---

def _custom_chunker(docs: List[Document], keyword_phrase: str = "List of all available keywords") -> List[Document]:
    """
    Chunks documents with a special rule for a keyword phrase.
    If a document contains the phrase, a single chunk is created from that phrase
    to the end of the document. The content before the phrase is chunked normally.
    """
    final_chunks = []
    standard_splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200)

    for doc in docs:
        content = doc.page_content
        metadata = doc.metadata

        start_index = content.find(keyword_phrase)

        if start_index != -1:
            # The phrase was found, apply the special rule

            # Part 1: Content before the phrase is split normally
            before_content = content[:start_index]
            if before_content.strip():
                before_chunks = standard_splitter.create_documents([before_content], metadatas=[metadata])
                final_chunks.extend(before_chunks)

            # Part 2: Content from the phrase to the end becomes one single chunk
            keyword_chunk_content = content[start_index:]
            keyword_chunk = Document(page_content=keyword_chunk_content, metadata=metadata)
            final_chunks.append(keyword_chunk)
        else:
            # Phrase not found, split the entire document normally
            chunks = standard_splitter.split_documents([doc])
            final_chunks.extend(chunks)

    return final_chunks

def load_all_docs_from_folder(folder_path: str,
         excluded_dirs: list[str] = None) -> (List[Document], List[str]): # Note the return type hint change
    """
    Loads all supported documents (.txt, .pdf, .html) from a directory,
    reliably excluding specified subdirectories and all their content.
    Returns a tuple of (all_docs, processed_files).
    """
    print(f"Loading documents from {folder_path}...")

    if excluded_dirs is None:
        excluded_dirs = set()
    else:
        excluded_dirs = set(excluded_dirs)

    all_docs = []
    processed_files = [] # Initialize a list for file paths ---
    loaders = {
        '.txt': TextLoader,
        '.pdf': PyPDFLoader,
        '.html': UnstructuredHTMLLoader
    }

    for dirpath, dirnames, filenames in os.walk(folder_path, topdown=True):
        dirnames[:] = [d for d in dirnames if d not in excluded_dirs]
        for filename in filenames:

            # Skip all hidden files (like .DS_Store)
            if filename.startswith('.'):
                continue

            file_path = os.path.join(dirpath, filename)
            file_ext = os.path.splitext(filename)[1].lower()

            if file_ext in loaders:
                loader_cls = loaders[file_ext]
                try:
                    loader = loader_cls(file_path)
                    all_docs.extend(loader.load())
                    processed_files.append(file_path) # --- Record the path
                    print("LOADING: %s (%s)" %(file_path, len(all_docs)))
                except Exception as e:
                    print(f"Error loading file {file_path}: {e}")

    print(f"Loaded {len(all_docs)} documents.")
    return all_docs, processed_files

def load_specific_docs(file_path_list: List[str]) -> List[Document]:
    """
    Loads a specific list of documents from their file paths.

    Args:
        file_path_list: A list of full paths to the files to be loaded.

    Returns:
        A list of LangChain Document objects.
    """
    all_docs = []
    loaders = {
        '.txt': TextLoader,
        '.pdf': PyPDFLoader,
        '.html': UnstructuredHTMLLoader
    }

    print(f"Loading {len(file_path_list)} specific files...")

    for file_path in file_path_list:
        # Determine the file extension
        file_ext = os.path.splitext(file_path)[1].lower()

        # If we have a loader for this extension, use it
        if file_ext in loaders:
            loader_cls = loaders[file_ext]
            try:
                # Initialize the loader with the specific file path
                loader = loader_cls(file_path)
                # Load the document and add its content to our list
                all_docs.extend(loader.load())
                print(f"  - Successfully loaded: {os.path.basename(file_path)}")
            except Exception as e:
                # Print an error if a specific file fails to load but continue
                print(f"  - WARNING: Error loading file {file_path}: {e}")
        else:
            print(f"  - WARNING: No loader available for file type: {file_path}")

    print(f"Successfully loaded content from {len(all_docs)} documents.")
    return all_docs

def get_llm_and_embeddings(
    provider: str = "google",
    llm_model_name: str = None,
    embedding_model_name: str = None,
    temperature: float = 0.1,
    timeout: int = 60,
    batch_size: int = 100
):
    """
    Initializes and returns the LLM and Embeddings clients for a given provider.

    Args:
        provider (str): "google" or "openai". Defaults to "google".
    """
    provider = provider.lower()
    if provider == "google":
        if llm_model_name is None:
            llm_model_name = "gemini-2.5-flash"
        if embedding_model_name is None:
            embedding_model_name = "models/embedding-001"

        llm = ChatGoogleGenerativeAI(
            model=llm_model_name,
            temperature=temperature,
            timeout=timeout,
            max_retries=0
        )
        embeddings = GoogleGenerativeAIEmbeddings(
            model=embedding_model_name,
            timeout=timeout,
            batch_size=batch_size,
            google_api_key=os.getenv("GOOGLE_API_KEY")
        )
        print("Using Google Gemini models.")

    elif provider == "openai":
        if llm_model_name is None:
            # Use a powerful and cost-effective model like gpt-4o
            llm_model_name = "gpt-4o"
        if embedding_model_name is None:
            # Use OpenAI's latest embedding model
            embedding_model_name = "text-embedding-3-small"

        llm = ChatOpenAI(
            model=llm_model_name,
            temperature=temperature,
            timeout=timeout,
            max_retries=2 # OpenAI client handles retries well
        )
        embeddings = OpenAIEmbeddings(
            model=embedding_model_name,
            chunk_size=batch_size # Corresponds to batch_size
        )
        print("Using OpenAI models.")

    else:
        raise ValueError("Unsupported provider. Choose 'google' or 'openai'.")

    return llm, embeddings

async def get_log_info(text, llm, embeddings, timeout: int = 120,
    provider: str = 'google'):
    """
    Stage 1 of the two-stage log file analysis pipeline with reranking.
    This stage summarizes the text from a log file and returns a dict with
    summary information.
    """
    try:
        # Pass the timeout to the summarization function
        log_summary_info = await summarize_log_text(
            text, llm, timeout=timeout, provider = provider)
        # process the log file to get pieces of information
        processed_log_dict = get_processed_log_dict(
           log_summary_info.log_summary)

        if log_summary_info.error:  # failed
           return group_args(
               group_args_type = 'error', error=log_summary_info.error)

        return group_args(group_args_type = 'log summary',
          summary = log_summary_info.log_summary,
          summary_as_html = save_as_html(log_summary_info.log_summary),
          processed_log_dict = processed_log_dict,
          error = None)

    except asyncio.TimeoutError:
        error_message = "Analysis timed out, " + \
             "try increasing timeout in Preferences or AnalyzeLog "+\
               "(currently %s sec)." %(timeout)

        print(error_message)
        return group_args(group_args_type = 'error', error=error_message)

    # --- OpenAI Specific Errors ---
    except openai.AuthenticationError:
        error_message = (
                "OPENAI API key is invalid"
            )

    except Exception as e:
        # Inspect the exception's original cause to find the specific error
        original_cause = e.__cause__

        if isinstance(original_cause, google_exceptions.PermissionDenied):
            error_message = (
                "Google AI API key is invalid"
            )
        elif isinstance(original_cause, google_exceptions.ResourceExhausted):
            error_message = (
                "Google AI API quota exceeded "
            )
        else:
            # Handle any other general exception
            if (str(e).find("API_KEY_INVALID")> -1) or (
                str(e).find("API_KEY_IP_ADDRESS_BLOCKED") > -1):
               error_message = "Google API key is invalid "

            elif str(e).find("free_tier_requests, limit: 0")> -1:
               error_message = "Google API key is not activated"

            elif str(e).find("limit: 0") > -1:
                  error_message = ( "Google AI API key has a zero quota" )



            else:
              error_message = f"An unexpected error occurred during summarization: {e}"
            print("Summarize log failed")

        print(error_message)
        return group_args(group_args_type='error', error=error_message)

async def analyze_log_summary(log_info, llm, embeddings,
   db_dir: str = "./docs_db",
   timeout: int = 60):
  # --- Stage 2: Analyze with Docs RAG

  try:
    vectorstore = load_persistent_db(embeddings, db_dir = db_dir)
    analysis_prompt = get_log_analysis_prompt()
    retriever = create_reranking_retriever(vectorstore, llm,
      timeout = timeout)
    analysis_rag_chain = create_log_analysis_chain(retriever, llm,
        analysis_prompt)
    retriever_query = ("Here is a summary of the %s log file:\n\n " %(
       log_info.processed_log_dict['phenix_program']) +
       log_info.processed_log_dict['summary'] +
       "\n\nConsidering whether the input data are from crystallography "+
       "(X-ray or neutron) or from cryo-EM, and considering the "+
       "the normal procedure for structure determination "
       "in Phenix, what are the next steps that I should carry out?" +
       "Consider this question in the context of the process of structure "+
       "determination in Phenix. Focus on using Phenix tools, but include "+
       "the use of Coot or Isolde if appropriate. Name the tools that are "+
       "to be used, along with their inputs and outputs and what they do."
      )

    final_analysis = analysis_rag_chain.invoke({
        "input": retriever_query,
        "log_summary": log_info.summary,
    })


    return group_args(group_args_type = 'answer',
      analysis = final_analysis.content,
      error = None)

# --- EXCEPTION HANDLERS ---
  except (TimeoutError, google_exceptions.DeadlineExceeded) as e:
      error_message = (
          "Network timeout with Google,"
          "You might try increasing the timeout in Preferences or in "
          f"AnalyzeLog (currently {timeout} sec)."
      )
      print(error_message)
      print(e)
      return group_args(group_args_type = 'answer',
        analysis = None,
        error = error_message)

  except google_exceptions.ResourceExhausted as e:
      error_message = (
          "ERROR: Google AI API quota exceeded."
          f"Details: {e}"
      )
      print(error_message)
      return group_args(group_args_type = 'answer',
        analysis = None,
        error = error_message)

  except CohereApiError as e:
      # Check if the error is specifically an authentication error (401)
      if hasattr(e, 'http_status') and e.http_status == 401:
            error_message = (
            "Invalid Cohere API key.\n"
            )
      elif str(e).find("invalid api token") > -1:
          error_message = "Cohere API key is invalid"
      else:
            # For other Cohere errors (like rate limiting)
            error_message = f"ERROR: A Cohere API error occurred, This may be due to rate limits.\nDetails: {e}"

      print(error_message)
      return group_args(group_args_type = 'answer',
        analysis = None,
        error = error_message)

  except Exception as e:
      error_message = "Reranking failed..."
      print(error_message)
      return group_args(group_args_type = 'answer',
        analysis = None,
        error = error_message)


def _iter_batches(seq: List[Document], size: int) -> Iterable[List[Document]]:
    for i in range(0, len(seq), size):
        yield seq[i:i+size]

def create_and_persist_db(
    docs: List[Document],
    embeddings,                       # GoogleGenerativeAIEmbeddings
    db_dir: str = "./docs_db",
    *,
    add_batch_size: int = 200,        # docs per add() call
    pause_between_batches: float = 2.0,
    max_attempts: int = 6,            # retries on rate limits
    max_backoff: float = 60.0
) -> Chroma:
    """
    Build a Chroma vector store with batching/throttling and automatic persistence
    (via chromadb.PersistentClient). No vectorstore.persist() needed in Chroma ≥0.5.
    """

    # Your existing chunking
    docs_chunks = _custom_chunker(docs)

    # Optional: light dedup to reduce embedding calls
    seen = set()
    unique_docs: List[Document] = []
    for d in docs_chunks:
        key = (d.page_content.strip(), tuple(sorted((d.metadata or {}).items())))
        if key not in seen:
            seen.add(key)
            unique_docs.append(d)

    # IMPORTANT: Use PersistentClient for on-disk storage
    client = chromadb.PersistentClient(path=db_dir)

    # Initialize an (initially empty) collection-backed vectorstore
    vectorstore = Chroma(
        client=client,
        collection_name="docs",
        embedding_function=embeddings,   # NOTE: embedding_function (not "embedding")
    )

    # Add in batches with backoff on 429/rate-limit errors
    for batch in _iter_batches(unique_docs, add_batch_size):
        backoff = 2.0
        for attempt in range(1, max_attempts + 1):
            try:
                vectorstore.add_documents(batch)
                break
            except Exception as e:
                msg = str(e).lower()
                if "API_KEY_IP_ADDRESS_BLOCKED" in msg:
                     error_message = (
                      "Google AI API key does not allow access from this server"
                       )
                     print(error_message)

                elif "you exceeded your current quota" in msg:
                   if "limit: 0" in msg:
                     error_message = (
                       "Google AI API key has a zero quota, "
                       )
                     print(error_message)
                     return None
                   else:
                     error_message = (
                       "Google AI API quota exceeded, "
                       )
                     print(error_message)
                     return None


                elif "429" in msg or "rate" in msg:
                    time.sleep(backoff)
                    backoff = min(backoff * 2.0, max_backoff)
                    if attempt == max_attempts:
                        raise RuntimeError("Failed to use Google API")
                else:
                    raise RuntimeError("Failed to use Google API")
        time.sleep(pause_between_batches)

    print(f"Database stored in: {db_dir} (persistence handled by PersistentClient)")
    return vectorstore

def load_persistent_db(embeddings: GoogleGenerativeAIEmbeddings,
      db_dir: str = "./docs_db") -> Chroma:
    """Loads a persisted Chroma vector store from disk."""
    if not os.path.exists(db_dir):
        raise FileNotFoundError(f"Database directory not found at '{db_dir}'.")

    # Load Chroma vector store from disk
    return Chroma(persist_directory=db_dir, embedding_function=embeddings)

def create_reranking_retriever(vectorstore: Chroma,
         llm: ChatGoogleGenerativeAI,
         timeout: int = 60):
    """
    Creates and returns a powerful retriever that uses a Cohere Reranker.
    """

    # Create the base retriever to get a large number of initial results
    base_retriever = vectorstore.as_retriever(search_kwargs={"k": 20})

    # 1. First, create a configured Cohere client with the API key and timeout.

    # Use the ClientV2 class as specified in the error message.
    cohere_client = cohere.ClientV2(
        api_key=os.getenv("COHERE_API_KEY"),
        timeout=timeout
    )

    # 2. Then, pass the pre-configured client to the CohereRerank object.
    reranker = CohereRerank(client=cohere_client, model="rerank-english-v3.0")
    #

    # The compression retriever combines the base retriever and the reranker
    compression_retriever = ContextualCompressionRetriever(
        base_compressor=reranker, base_retriever=base_retriever
    )

    return compression_retriever

def create_reranking_rag_chain(retriever, llm: ChatGoogleGenerativeAI,
        prompt: PromptTemplate):
    """
    Creates a full RAG chain by combining a reranking retriever and a prompt.
    """
    # Set up RAG chain with Cohere Reranker

    def format_docs(docs):
        return "\n\n".join(doc.page_content for doc in docs)

    # Build the final LCEL chain
    rag_chain = (
        {"context": retriever | format_docs, "input": RunnablePassthrough()}
        | prompt
        | llm
    )
    return rag_chain

def create_log_analysis_chain(retriever, llm: ChatGoogleGenerativeAI, prompt: PromptTemplate):
    """
    Creates the specific RAG chain for log analysis.
    """
    # Set up log analysis RAG chain

    # Define the inputs for the final prompt
    inputs = {
        "context": lambda x: retriever.invoke(x["input"]),
        "log_summary": lambda x: x["log_summary"]
    }

    # Build the final chain
    analysis_rag_chain = (
        RunnablePassthrough.assign(**inputs)
        | prompt
        | llm
    )
    return analysis_rag_chain

def get_processed_log_dict(log_text: str ) -> dict:
  # Break up log_text into sections based on the text like:
  #  **8. Key Output Files and Evaluation:**
  summary = find_text_block(log_text, "**8", end_text = "**")
  phenix_program = find_text_block(log_text, "**3", end_text = "**")

  dd = {'phenix_program':phenix_program,
        'summary':summary}
  return dd

def find_text_block(log_text: str, target_text: str, end_text: str = "**"):
  text_lines = []
  started = False
  for line in log_text.splitlines():
    if (not started) and line.strip().startswith(target_text):
       started = True
    elif (started) and line.strip().startswith(end_text):
       break
    if started:
      text_lines.append(line)
  return " ".join(text_lines)

def query_docs(query_text, llm=None, embeddings=None,
        db_dir: str = "./docs_db",
        timeout: int = 60,
        max_attempts: int = 5, # Maximum number of retries
        use_throttling: bool = False,
        provider: str = "google"):
    """Query Phenix docs with query_text,
    with automatic retries on rate limits."""

    global _last_query_time

    # --- Conditional Throttling Block ---
    if use_throttling:
        min_interval = 4  # seconds, for ~15 RPM
        elapsed_time = time.time() - _last_query_time
        if elapsed_time < min_interval:
            wait_time = min_interval - elapsed_time
            print(f"Throttling: Waiting {wait_time:.1f} seconds ...")
            time.sleep(wait_time)

        # Update the timestamp only when throttling is active
        _last_query_time = time.time()
    # ------------------------------------------

    if not llm:
      # Set up the LLM and embeddings
      try:
        llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout)
      except ValueError as e:
        print(e)
        raise ValueError("Sorry, unable to set up LLM with %s" %(provider))



    # Add some qualifications to the query text to make it focus on Phenix
    query_text += """Consider this question in the context of the process
        of structure determinaion in Phenix. Focus on using Phenix tools,
        but include the use of Coot or Isolde if appropriate. Name the tools
        that are to be used, along with their inputs and outputs and what
        they do."""

    # Get the data from the database
    vectorstore = load_persistent_db(embeddings, db_dir=db_dir)

    # Set up retriever
    retriever = create_reranking_retriever(vectorstore, llm, timeout=timeout)

    # Get the prompt
    prompt = get_docs_query_prompt()

    #  Create the advanced chain
    rag_chain = create_reranking_rag_chain(retriever, llm, prompt)

    # --- NEW: Add retry loop with exponential backoff ---
    backoff_time = 2.0  # Initial wait time in seconds
    for attempt in range(max_attempts):
        try:
            # Get the response and return its content
            response = rag_chain.invoke(query_text)
            return response.content

        # --- OpenAI Specific Errors ---
        except openai.RateLimitError:
            if attempt < max_attempts - 1:
                print(f"OpenAI rate limit exceeded. Waiting {backoff_time:.1f} sec...")
                time.sleep(backoff_time)
                backoff_time *= 2
            else:
                print("OpenAI API rate limit exceeded after multiple retries.")
                return None
        except openai.AuthenticationError:
            print("OpenAI API key is invalid")
            return None


        except (google_exceptions.ResourceExhausted,
                 GoogleGenerativeAIError) as e:
            # check the error message text
            error_text = str(e).lower()
            if "you exceeded your current quota" in error_text:
                if "limit: 0" in error_text:
                  error_message = (
                    "Google AI API has a zero quota, "
                    "Please go to console.cloud.google.com if this is"
                    "your plan and check if it is active.\n"
                    )
                  print(error_message)
                  return None
                else:
                  error_message = (
                    "Google AI API quota exceeded, "
                    "Please go to console.cloud.google.com if this is"
                    "your plan and check if it is active.\n"
                    )
                  print(error_message)
                  return None

            if "429" in error_text or "rate limit" in error_text:
                if attempt < max_attempts - 1:
                    print("Rate limit exceeded. ")
                    print(f"Waiting for {backoff_time:.1f} sec...")
                    time.sleep(backoff_time)
                    backoff_time *= 2
                else:

                  error_message = (
                    "Google AI API quota exceeded, "
                    "Please check your plan and billing details.\n"
                    )
                  print(error_message)
                  return None
            else:
                # If it's a different kind of GoogleGenerativeAIError
                #    (not a rate limit), handle it here
                print(f"An unexpected Google AI error occurred: {e}")
                return None


        except Exception as e:
            # Handle other exceptions
            print(f"An unexpected error occurred during query: {e}")
            return None

    # This part is reached if all retries fail
    print("Query docs failed...")
    return None


def get_chunk_size(provider: str = 'google'):
    provider = provider.lower()
    if provider == "openai":
        chunk_size = 100000  # Safe for GPT-4o's 128k limit
        chunk_overlap = 10000
        print("Using smaller chunk size for OpenAI.")
    elif provider == "google":
        chunk_size = 750000  # Safe for Gemini 2.5 Flash 1M limit
        chunk_overlap = 50000
        print("Using larger chunk size for Google Gemini.")
    else:
        # Default to the safest (smallest) size if provider is unknown
        chunk_size = 100000
        chunk_overlap = 10000
        print(f"Warning: Unknown provider '{provider}'. Defaulting to smaller chunk size.")
    # -------------------------------------
    return chunk_size, chunk_overlap

def _custom_log_chunker(log_text: str,
     provider: str = "google") -> List[Document]:
    """
    Chunks a log text with a special rule for the 'Files are in the directory' section.

    - If the trigger phrase is found, the text before it is chunked normally.
    - The text from the trigger phrase to the end of the file (or to 'Citations')
      is treated as a single, separate chunk.
    - If the trigger phrase is not found, the entire log is chunked normally.
    """


    chunk_size, chunk_overlap = get_chunk_size(provider)

    trigger_phrase = "Files are in the directory"
    end_phrase = "Citations"

    final_chunks = []

    # Search for the trigger phrase case-insensitively
    start_index = log_text.lower().find(trigger_phrase.lower())

    if start_index == -1:
        # Trigger not found, use the standard chunker for the whole document
        documents = [Document(page_content=log_text)]
        standard_splitter = RecursiveCharacterTextSplitter(
            chunk_size=chunk_size, chunk_overlap=chunk_overlap
        )
        return standard_splitter.split_documents(documents)

    # --- Trigger was found, apply special logic ---

    # 1. Chunk the text *before* the trigger phrase normally
    before_text = log_text[:start_index].strip()
    if before_text:
        before_doc = [Document(page_content=before_text)]
        standard_splitter = RecursiveCharacterTextSplitter(
            chunk_size=chunk_size, chunk_overlap=chunk_overlap
        )
        final_chunks.extend(standard_splitter.split_documents(before_doc))

    # 2. Create the single, special chunk
    special_section_text = log_text[start_index:]

    # Find the end phrase within the special section
    end_index = special_section_text.find(end_phrase)

    if end_index != -1:
        # If 'Citations' is found, end the chunk there
        special_chunk_content = special_section_text[:end_index].strip()
    else:
        # Otherwise, the chunk goes to the end of the file
        special_chunk_content = special_section_text.strip()

    if special_chunk_content:
        final_chunks.append(Document(page_content=special_chunk_content))

    return final_chunks

async def summarize_log_text(
    text: str,
    llm: ChatGoogleGenerativeAI,
    timeout: int = 120,
    batch_size: int = 10,  # Process 10 chunks at a time (safely under 15 RPM limit)
    pause_between_batches: int = 1,  # Wait 1 seconds between batches
    use_throttling: bool = False,
    provider: str = 'google',
) -> str:
    """
    Performs a map-reduce summarization with batching to respect API rate limits.
    """
    docs = _custom_log_chunker(text, provider = provider)
    if not docs:
      return group_args(group_args_type = 'log_summary',
        log_summary = None,
        error = "Log file produced no content to summarize.")

    map_prompt = get_log_map_prompt()
    map_chain = map_prompt | llm

    all_intermediate_summaries = []

    # Use the existing _iter_batches helper to process chunks in batches
    num_batches = (len(docs) + batch_size - 1) // batch_size
    print(f"Summarizing {len(docs)} chunks in {num_batches} batches to respect rate limits...")

    for i, batch in enumerate(_iter_batches(docs, batch_size)):
        print(f"  - Processing batch {i + 1} of {num_batches}...")

        # Create tasks for the current batch only
        tasks = [map_chain.ainvoke({"text": doc.page_content}) for doc in batch]

        # Await the completion of the current batch
        try:
            map_results = await asyncio.wait_for(asyncio.gather(*tasks), timeout=timeout)
            all_intermediate_summaries.extend([result.content for result in map_results])
        except asyncio.TimeoutError:
            # Handle timeout for the batch
            print(f"Batch {i + 1} timed out after {timeout} seconds. Proceeding with partial results.")
            continue # You could also choose to raise the error here

        # If it's not the last batch and throttling is enabled, pause.
        if use_throttling and (i < num_batches - 1):
            print(f"  - Pausing {pause_between_batches} seconds ...")
            await asyncio.sleep(pause_between_batches)
        # ----------------------------

    # --- The rest of the function proceeds as before with the collected summaries ---
    if not all_intermediate_summaries:
      return group_args(group_args_type = 'log_summary',
        log_summary = None,
        error = "Log summarization failed to produce any results.")

    summary_docs = [Document(page_content=s) for s in all_intermediate_summaries]

    combine_prompt = get_log_combine_prompt()
    reduce_chain = create_stuff_documents_chain(llm, combine_prompt)
    final_output = reduce_chain.invoke({"context": summary_docs})
    return group_args(group_args_type = 'log_summary',
      log_summary = final_output,
      error = None)
