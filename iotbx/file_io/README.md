# iotbx.file_io

Fast, robust **file-type identification** and a **parse-once read facade** for the
`iotbx.data_manager` datatypes. This package answers two questions about a file:

1. *What is it?* — `get_file_type(path)` returns a DataManager datatype
   (`model`, `miller_array`, `real_map`, `restraint`, `sequence`, `phil`,
   `ncs_spec`, `json`, `yaml`) **without a full parse**, or `None` if it can't
   be identified. It never raises.
2. *Load it* — `read_file(path)` parses the file **exactly once**, with the type
   already known, and returns a uniform `FileIOResult`.

`get_file_type` is the modern, Qt-independent replacement for using
`iotbx.file_reader.any_file` just to learn a type; `read_file` replaces using
`any_file(path).file_object` to load it.

## Public API

```python
from iotbx.file_io import get_file_type, read_file, FileIOResult

get_file_type(path, valid_types=None, verify=True, logger=None, cif_engine='xcif')
    # -> datatype str | None  (never raises)

read_file(path, file_type=None, cif_engine='xcif', force=False)
    # -> FileIOResult(data_type, file_object, reader)  (raises Sorry on failure)
```

- `valid_types` — restrict the answer to a set of datatypes; a detected type not
  in the set yields `None`. (The DataManager passes its own `datatypes` here.)
- `verify` — when `True` (default), confirm the extension's guess against the
  file's content; when `False`, trust the extension.

Also re-exported: the maps `extension_to_datatype`, `any_file_type`,
`data_manager_type`, and the `CIF_SENTINEL` marker.

---

## How a file is identified — `get_file_type`

Identification is **tiered**, cheapest first, and only does as much work as
needed:

1. **Normalize** — strip a SHELX format suffix, then a compression-aware
   `splitext` into `(base, extension, compress_ext)`.
2. **Extension → candidate** datatype via `extension_to_datatype`
   (`.pdb`→`model`, `.ncs`→`ncs_spec`, `.eff`→`phil`, `.cif`→*CIF marker*, …).
3. **Compression gate** — decide whether a cheap decompressed prefix peek is
   possible (`.gz/.bz2/.xz/.lzma/.zst/.zstd` yes; `.Z` no; `.zip` is unreadable).
4. **CIF** (`.cif`/`.mmcif`) — resolved by an authoritative but fast `xcif`
   token parse into `model` / `miller_array` / `restraint` (and combined CIFs).
5. **Cheap content verification** of a bounded (≤64 KB) decompressed prefix:
   - **binary magic** for `.mtz` (`MTZ `) and CCP4/MRC maps (`MAP ` at byte 208);
   - **text-vs-binary** gate, then a **structural content sniff** for the text
     types that can collide (see *Reclassifying* below).
6. **`any_file` fallback** — whenever the cheap tiers are inconclusive, defer to
   a full `any_file` parse, so detection is never *less* capable than `any_file`.
7. **Filter** by `valid_types`.

The whole body is wrapped so the **never-raises** contract holds: expected
errors (missing/unreadable file, decompression failure, parse error) become
`None`; a truly unexpected error is logged (or re-raised when the
`IOTBX_FILE_IO_DEBUG` environment variable is set) and then degrades to `None`.

### Flow

```
 get_file_type(path, valid_types=None, verify=True)
   |
   |-- strip SHELX suffix; splitext -> (base, ext, compress_ext)
   |
   |-- path is not a file? --------------------------------------> None   ........ FAIL  (missing/unreadable)
   |
   |-- candidate = extension_to_datatype[ext]      # .pdb->model, .ncs->ncs_spec, ...
   |
   |-- compression gate (compress_ext):
   |       .zip ------------------------------------------------> None   ........ FAIL  (not readable)
   |       .Z  / unknown compressor ......................... can_peek = False
   |       .gz .bz2 .xz .lzma .zst .zstd .................... can_peek = True
   |
   |== candidate == CIF_SENTINEL   (.cif / .mmcif)
   |       xcif token-parse -> model | miller_array | restraint -> datatype  .... WORK   (CIF by content)
   |       inconclusive ---------------------------------------> any_file fallback
   |
   |== candidate is None   (unknown extension) ----------------> any_file fallback
   |
   |== not verify  OR  not can_peek   (.Z, missing decompressor)
   |       -----------------------------------------------------> _filtered(candidate)  WORK  (extension-trust)
   |
   |-- read <=64 KB decompressed prefix
   |
   |== BINARY MAGIC
   |       .mtz  : leading "MTZ "  ?   yes -> miller_array  .... WORK     no -> any_file fallback
   |       real_map : "MAP " @ byte 208 ?  yes -> real_map  ... WORK     no -> any_file fallback
   |
   |== TEXT GATE   (candidate in model/sequence/phil/ncs_spec/restraint/json/yaml)
   |       content is binary?  --------------------------------> any_file fallback  ... (misnamed binary)
   |       |
   |       sniffable text (model / sequence / phil / ncs_spec):
   |          sniff_text_datatype(prefix):
   |             match    -> _filtered(SNIFFED) <=============== RECLASSIFY  (content wins over extension)
   |             no match -> any_file fallback
   |       |
   |       json / yaml (no sniffer) ----------------------------> _filtered(candidate)  WORK  (extension-trust)
   |
   |== otherwise (hkl/sca/... non-mtz miller_array) -----------> _filtered(candidate)  WORK  (extension-trust)


 any_file fallback : data_manager_type[ any_file(path).file_type ]   (full parse)  -> datatype | None
 _filtered(t)      : t  if (valid_types is None or t in valid_types)  else None   ... FAIL  (valid_types)
 error wrapper     : I/O / decompression / parse error -> None ;
                     unexpected error -> warn (or raise if IOTBX_FILE_IO_DEBUG) -> None
```

---

## The three kinds of path

### 1. Working (happy) path

The extension names a datatype and the content agrees, so identification stays on
the cheap path — no full parse:

- **By magic bytes** — `model.mtz` (`MTZ ` stamp) → `miller_array`;
  `map.ccp4` (`MAP ` at offset 208) → `real_map`.
- **By structural sniff** — `model.pdb` (an `ATOM`/`HETATM` record) → `model`;
  `seq.fasta` (residue letters) → `sequence`; `params.eff`
  (`name = value` / `scope {`) → `phil`; `ops.ncs_spec` (`new_ncs_group`,
  `rota_matrix`, …) → `ncs_spec`.
- **By CIF token-scan** — `lig.cif` → `restraint`, `model.cif` → `model`,
  `data.cif` → `miller_array`, decided from the categories present.
- **By extension-trust** — `json`/`yaml` (which `any_file` can't classify),
  the non-`mtz` reflection extensions (`.hkl`, `.sca`, …), and any file read with
  `verify=False` or a `.Z`/non-peekable compressor are accepted on their
  extension alone.

Compression is transparent: the prefix is decompressed before any check, so
`model.pdb.gz`, `lig.cif.xz`, etc. identify exactly like their uncompressed
forms.

### 2. Reclassifying path — *content wins over the extension*

The four **sniffable text types** — `model` (PDB), `sequence`, `phil`,
`ncs_spec` — can be mislabeled by extension because their formats overlap or
their extensions are weak. After the text-vs-binary gate passes,
`sniff_text_datatype` reads the prefix and returns the type the **content**
actually matches, in precedence order `model → ncs_spec → phil → sequence`:

- A **PHIL** refinement override saved as `override.ncs` → detected `phil`
  (not `ncs_spec`). This is the motivating case: without it the file is taken for
  an `ncs_spec`, fails to load, and is silently dropped.
- A PHIL body misnamed `override.pdb` → `phil`.
- The sniffers are tuned to agree with what the format reader will actually
  accept (e.g. PDB requires an atom-bearing record; a *gapped* `.seq` alignment
  is **not** `sequence`), so a confident sniff never disagrees with a successful
  later parse.

If no sniffer matches, detection defers to the `any_file` fallback rather than
trusting the extension. `json`/`yaml` are **not** reclassified — they have no
sniffer and keep extension-trust (`any_file` cannot classify them).

### 3. Failure paths

`get_file_type` returns `None` (never raises) when:

| Situation | Outcome |
|---|---|
| File missing / unreadable | `None` |
| `.zip` archive | `None` (not readable by the I/O layer) |
| Corrupt or **truncated** compressed file | `None`, silently — decompression errors (`EOFError`, `zlib.error`, `lzma.LZMAError`, `zstandard.ZstdError`) are caught, matching a corrupt `.gz` |
| Unknown extension **and** unrecognized content | `None` (via the `any_file` fallback) |
| Misnamed binary under a text extension | re-routed by the `any_file` fallback → correct type or `None` |
| Detected type excluded by `valid_types` | `None` |
| Unexpected internal error | `None`, after a warning (or re-raised with `IOTBX_FILE_IO_DEBUG`) |

At the **load** layer, `read_file` does raise `Sorry` for a missing file, an
unrecognized type, or a parse failure of the chosen type.

---

## Reading a file — `read_file`

`read_file` is the single place a data file is parsed for content. With the type
known (passed in, or detected via `get_file_type`) it dispatches once to the
existing format reader and wraps the result:

| datatype | reader | `FileIOResult.file_object` |
|---|---|---|
| `model` | `iotbx.pdb.hierarchy.input` (PDB + mmCIF) | pdb hierarchy input (`.input` available) |
| `miller_array` / `map_coefficients` | `iotbx.reflection_file_reader.any_reflection_file` | reflection file object |
| `real_map` | `iotbx.map_manager` | map object |
| `sequence` | `iotbx.bioinformatics.any_sequence_format` | sequence objects |
| `ncs_spec` | `mmtbx.ncs` (BIOMT records, else the ncs_spec text reader) | ncs object |
| `phil` | `iotbx.phil.parse` | phil object |
| `restraint` | `iotbx.cif.reader(engine=cif_engine)` | cif reader (`.model()`) |
| `json` / `yaml` | `json.load` / `yaml.safe_load` | parsed Python object |

The datatypes that were once `any_file`-backed now read **directly** via the
readers above; `any_file` is retained only as a fallback during the deprecation
window (a direct read leaves `FileIOResult.reader is None`, while the fallback
sets it to the `any_file` input). When a reader auto-detects, `read_file` accepts
the parse only if the detected type matches (raising `Sorry` on a mismatch);
`force=True` parses as the requested type instead — the mechanism that pulls a
single datatype out of a combined CIF.

A compressed file whose underlying reader cannot read the compressor directly is
transparently decompressed to a temp file and retried; `.Z` is decompressed via
`gunzip`. `FileIOResult.file_content` is an alias for `file_object`
(mirroring `any_file`).

`DataManager.process_file(path)` ties the two layers together: it calls
`get_file_type`, then the matching `process_<datatype>_file` (which uses
`read_file`). If detection picks a type but the reader rejects the content, it
reports the file as unused (`[]`) rather than raising.

---

## Module layout

| File | Responsibility |
|---|---|
| `__init__.py` | public API: `get_file_type`, `read_file`, `FileIOResult`, the maps |
| `detection.py` | `get_file_type` — the tiered identification procedure above |
| `text_detection.py` | `sniff_text_datatype` + the structural text sniffers (`_looks_like_pdb/_ncs_spec/_phil/_sequence`); the reclassifying engine |
| `reader.py` | `read_file` (the single parse authority) + `FileIOResult` |
| `_maps.py` | datatype ↔ `any_file` maps + the extension→datatype table |

## Extension map (Tier 1)

| extension(s) | candidate datatype |
|---|---|
| `pdb`, `ent` | `model` |
| `mtz`, `hkl`, `sca`, `cns`, `cv`, `ref`, `fobs` | `miller_array` |
| `cif`, `mmcif` | *CIF marker* → resolved by content |
| `ccp4`, `mrc`, `map` | `real_map` |
| `ncs`, `ncs_spec` | `ncs_spec` |
| `params`, `eff`, `def`, `phil`, `param` | `phil` |
| `fa`, `faa`, `seq`, `pir`, `dat`, `fasta` | `sequence` |
| `json` | `json` |
| `yaml`, `yml` | `yaml` |

Compression suffixes (`.gz`, `.bz2`, `.xz`, `.lzma`, `.zst`, `.zstd`, `.Z`) are
stripped by the compression-aware `splitext` before this lookup.
