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
import cohere
from cohere.core.api_error import ApiError as CohereApiError

import os
os.environ['GRPC_ENABLE_FORK_SUPPORT'] = "false"   # suppress a warning message

import asyncio
from concurrent.futures import TimeoutError
from libtbx import group_args
from google.api_core import exceptions as google_exceptions
from langchain_cohere import CohereRerank
from langchain.retrievers import ContextualCompressionRetriever
from langchain_core.runnables import RunnablePassthrough
from langchain_google_genai import ChatGoogleGenerativeAI, GoogleGenerativeAIEmbeddings
from langchain_community.document_loaders import TextLoader, PyPDFLoader, UnstructuredHTMLLoader
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.chains.combine_documents import create_stuff_documents_chain
from langchain_core.prompts import PromptTemplate
from langchain_core.documents import Document
from langchain_chroma import Chroma
from typing import List
from markdown_it import MarkdownIt

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
    4.  **Be Concise:** Create a brief, bulleted list of these key facts. Do not add conversational text or explanations.

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
        "1. Table of input files and their contents.\n"
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
        "Based on all the information above, please perform the following analysis:\n"
        "1. Place the events of the log summary in the broader context of a "
        "typical Phenix structure determination workflow as described in "
        "the documentation and papers."
        "For each step in the workflow, mention the inputs and outputs. "
        "If any metrics are available for that step, specifically mention "
        "those metrics and comment on the quality of the result based on "
        "those metrics. "

        "2. List the key output files, along with the values of any available "
         "metrics describing their utilities.  If no metrics are available, "
         "do not provide any. "
        "3. Evaluate whether the run described in the summary was useful. "
          "List reported metrics and expected values of these metrics and "
          "consider the goals of the program.\n"
        "4. Considering the normal sequence of Phenix tool use, suggest "
        "three concrete next steps in structure determination using Phenix "
        "or graphical tools such as Coot or Isolde that I should take, "
        "justifying each suggestion with information from the provided "
        "documentation and papers. "
        "Do not include any information on CryoFit, ShelxD, or CCP4 tools "
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
         excluded_dirs: list[str] = None) -> List[Document]:
    """Loads all supported documents (.txt, .pdf, .html) from a directory,
        excluding specified subdirectories."""
    print(f"Loading documents from {folder_path}...")

    if excluded_dirs:
      exclude_patterns = [f"{d}/**" for d in excluded_dirs]
    else:
      exclude_patterns = []

    loaders = {
        '.txt': TextLoader,
        '.pdf': PyPDFLoader,
        '.html': UnstructuredHTMLLoader
    }

    all_docs = []

    # os.walk is a standard Python function to recursively explore a directory.
    for dirpath, dirnames, filenames in os.walk(folder_path, topdown=True):

        # We modify 'dirnames' in-place. os.walk will see this modification
        # and will NOT recurse into any of the directories we remove.
        dirnames[:] = [d for d in dirnames if d not in excluded_dirs]

        for filename in filenames:
            # Construct the full path to the file
            file_path = os.path.join(dirpath, filename)

            # Determine the file extension
            file_ext = os.path.splitext(filename)[1].lower()

            # If we have a loader for this extension, use it
            if file_ext in loaders:
                loader_cls = loaders[file_ext]
                try:
                    # Initialize the loader with the specific file path
                    loader = loader_cls(file_path)
                    # Load the document and add its content to our list
                    all_docs.extend(loader.load())
                    print("LOADING: %s (%s)" %(file_path, len(all_docs)))
                except Exception as e:
                    # Print an error if a specific file fails to load
                    print(f"Error loading file {file_path}: {e}")


    print(f"Loaded {len(all_docs)} documents.")
    return all_docs

def get_llm(model_name: str = "gemini-1.5-pro-latest", temperature: float = 0.1, timeout: int = 60) -> ChatGoogleGenerativeAI:
    """Initializes and returns a ChatGoogleGenerativeAI instance."""
    return ChatGoogleGenerativeAI(model=model_name, temperature=temperature, timeout=timeout)

def get_embeddings(model_name: str = "models/embedding-001", timeout: int = 60) -> GoogleGenerativeAIEmbeddings:
    """Initializes and returns a GoogleGenerativeAIEmbeddings instance."""
    return GoogleGenerativeAIEmbeddings(model=model_name, timeout=timeout)

# In langchain_tools.py

async def get_log_info(text, llm, embeddings, timeout: int = 120):
    """
    Stage 1 of the two-stage log file analysis pipeline with reranking.
    This stage summarizes the text from a log file and returns a dict with
    summary information.
    """
    try:
        # Pass the timeout to the summarization function
        log_summary = await summarize_log_text(
            text, llm, timeout=timeout)
        # process the log file to get pieces of information
        processed_log_dict = get_processed_log_dict(log_summary)

        return group_args(group_args_type = 'log summary',
          summary = log_summary,
          summary_as_html = save_as_html(log_summary),
          processed_log_dict = processed_log_dict,
          error = None)

    except asyncio.TimeoutError:
        error_message = f"The Google analysis timed out.\n" +\
         "You might try increasing the timeout in Preferences or in " +\
          f"AnalyzeLog (currently {timeout} sec)"

        print(error_message)
        return group_args(group_args_type = 'error', error=error_message)

    except Exception as e:
        # Inspect the exception's original cause to find the specific error
        original_cause = e.__cause__

        if isinstance(original_cause, google_exceptions.PermissionDenied):
            error_message = (
                "Google AI Permission Denied. "
                 "This usually means your API key is invalid, "
                "has been disabled, or is not enabled for the "
                "Gemini API on your Google Cloud project.\n"
                f"Details: {original_cause}"
            )
        elif isinstance(original_cause, google_exceptions.ResourceExhausted):
            error_message = (
                "Google AI API quota exceeded. "
                "Please check your plan and billing details.\n"
                f"Details: {original_cause}"
            )
        else:
            # Handle any other general exception
            if str(e).find("API_KEY_INVALID")> -1:
               error_message = "The Google API key is invalid. "+ \
                 "Please check in Preferences"
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
       "\n\nConsidering the normal procedure for structure determination "
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
          "A network timeout occurred while communicating Google.\n"
          "You might try increasing the timeout in Preferences or in "
          f"AnalyzeLog (currently {timeout} sec)"
      )
      print(error_message)
      print(e)
      return group_args(group_args_type = 'answer',
        analysis = None,
        error = error_message)

  except google_exceptions.ResourceExhausted as e:
      error_message = (
          "ERROR: Google AI API quota exceeded."
          " Please check your plan and billing details.\n"
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
            "ERROR: Cohere API Authentication Error."
            " This almost always means your Cohere API key is "
            "invalid or missing.\n"
                f"Details: {e}"
            )
      elif str(e).find("invalid api token") > -1:
          error_message = "The Cohere API key is invalid. "+ \
                 "Please check in Preferences"
      else:
            # For other Cohere errors (like rate limiting)
            error_message = f"ERROR: A Cohere API error occurred. This may be due to rate limits.\nDetails: {e}"

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

def create_and_persist_db(docs: List[Document],
      embeddings: GoogleGenerativeAIEmbeddings,
      db_dir: str = "./docs_db") -> Chroma:
    """Creates and persists a Chroma vector store from documents."""

    docs_chunks = _custom_chunker(docs)

    # Create and save Chroma vector store to db
    vectorstore = Chroma.from_documents(
        documents=docs_chunks,
        embedding=embeddings,
        persist_directory=db_dir
    )
    # Load Chroma vector store from disk
    print("Database to be written to: %s" %(db_dir))
    return Chroma(persist_directory=db_dir, embedding_function=embeddings)


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

def query_docs(query_text, llm = None, embeddings = None,
        db_dir: str = "./docs_db",
        timeout: int = 60):
  """Query Phenix docs with query_text"""
  try:
    if not llm:
      # Set up the LLM and embeddings
      llm = get_llm()
      embeddings = get_embeddings()

    # Add some qualifications to the query text to make it focus on Phenix

    query_text += """Consider this question in the context of the process
        of structure determinaion in Phenix. Focus on using Phenix tools,
        but include the use of Coot or Isolde if appropriate. Name the tools
        that are to be used, along with their inputs and outputs and what
        they do."""


    # Get the data from the database
    vectorstore = load_persistent_db(embeddings, db_dir = db_dir)

    # Set up retriever
    retriever = create_reranking_retriever(vectorstore, llm,
      timeout = timeout)

    # Get the prompt
    prompt = get_docs_query_prompt()

    #  Create the advanced chain
    rag_chain = create_reranking_rag_chain(retriever, llm, prompt)

    # Get the response and return its content
    response = rag_chain.invoke(query_text)
    return response.content
  except Exception as e:
    print("Query docs failed...")
    return None

# Change the signature and add asyncio.wait_for
async def summarize_log_text(
    text: str, llm: ChatGoogleGenerativeAI, timeout: int = 120) -> str:
    """Performs the map-reduce summarization for a single log file."""
    # ... (code for creating documents and splitting text remains the same) ...
    documents = [Document(page_content=text)]
    text_splitter = RecursiveCharacterTextSplitter(chunk_size=300000,
         chunk_overlap=30000)
    docs = text_splitter.split_documents(documents)

    map_prompt = get_log_map_prompt()
    map_chain = map_prompt | llm

    tasks = [map_chain.ainvoke({"text": doc.page_content}) for doc in docs]

    # --- THIS IS THE KEY CHANGE ---
    # We wrap the asyncio.gather call in asyncio.wait_for
    map_results = await asyncio.wait_for(asyncio.gather(*tasks), timeout=timeout)
    # ----------------------------

    intermediate_summaries = [result.content for result in map_results]

    summary_docs = [Document(page_content=s) for s in intermediate_summaries]

    combine_prompt = get_log_combine_prompt()
    reduce_chain = create_stuff_documents_chain(llm, combine_prompt)
    final_output = reduce_chain.invoke({"context": summary_docs})
    return final_output

def save_as_html(markdown_string: str,
     title: str = "Summary", file_name: str = None):
    """Converts a Markdown string to an HTML file."""
    md = MarkdownIt("gfm-like")
    html_content = md.render(markdown_string)

    # Add some basic styling for better readability
    html_with_style = f"""
    <html>
    <head>
        <title>{title}</title>
        <style>
            body {{ font-family: sans-serif; line-height: 1.6; padding: 2em; max-width: 800px; margin: auto; }}
            table {{ border-collapse: collapse; width: 100%; margin-bottom: 1em; }}
            th, td {{ border: 1px solid #dddddd; text-align: left; padding: 8px; }}
            th {{ background-color: #f2f2f2; }}
            code {{ background-color: #eee; padding: 2px 4px; border-radius: 3px; }}
        </style>
    </head>
    <body>
        <h2></b>{title}</b></h2>
        {html_content}
    </body>
    </html>
    """

    if file_name:
      with open(file_name, "w") as f:
        f.write(html_with_style)
      print(f"\nSaved formatted output to: {file_name}")

    return html_with_style
