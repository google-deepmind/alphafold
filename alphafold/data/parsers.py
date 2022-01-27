# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Functions for parsing various file formats."""
import collections
import dataclasses
import itertools
import re
import string
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Set

# Internal import (7716).


DeletionMatrix = Sequence[Sequence[int]]


@dataclasses.dataclass(frozen=True)
class Msa:
  """Class representing a parsed MSA file."""
  sequences: Sequence[str]
  deletion_matrix: DeletionMatrix
  descriptions: Sequence[str]

  def __post_init__(self):
    if not (len(self.sequences) ==
            len(self.deletion_matrix) ==
            len(self.descriptions)):
      raise ValueError(
          'All fields for an MSA must have the same length. '
          f'Got {len(self.sequences)} sequences, '
          f'{len(self.deletion_matrix)} rows in the deletion matrix and '
          f'{len(self.descriptions)} descriptions.')

  def __len__(self):
    return len(self.sequences)

  def truncate(self, max_seqs: int):
    return Msa(sequences=self.sequences[:max_seqs],
               deletion_matrix=self.deletion_matrix[:max_seqs],
               descriptions=self.descriptions[:max_seqs])


@dataclasses.dataclass(frozen=True)
class TemplateHit:
  """Class representing a template hit."""
  index: int
  name: str
  aligned_cols: int
  sum_probs: Optional[float]
  query: str
  hit_sequence: str
  indices_query: List[int]
  indices_hit: List[int]


def parse_fasta(fasta_string: str) -> Tuple[Sequence[str], Sequence[str]]:
  """Parses FASTA string and returns list of strings with amino-acid sequences.

  Arguments:
    fasta_string: The string contents of a FASTA file.

  Returns:
    A tuple of two lists:
    * A list of sequences.
    * A list of sequence descriptions taken from the comment lines. In the
      same order as the sequences.
  """
  sequences = []
  descriptions = []
  index = -1
  for line in fasta_string.splitlines():
    line = line.strip()
    if line.startswith('>'):
      index += 1
      descriptions.append(line[1:])  # Remove the '>' at the beginning.
      sequences.append('')
      continue
    elif not line:
      continue  # Skip blank lines.
    sequences[index] += line

  return sequences, descriptions


def parse_stockholm(stockholm_string: str) -> Msa:
  """Parses sequences and deletion matrix from stockholm format alignment.

  Args:
    stockholm_string: The string contents of a stockholm file. The first
      sequence in the file should be the query sequence.

  Returns:
    A tuple of:
      * A list of sequences that have been aligned to the query. These
        might contain duplicates.
      * The deletion matrix for the alignment as a list of lists. The element
        at `deletion_matrix[i][j]` is the number of residues deleted from
        the aligned sequence i at residue position j.
      * The names of the targets matched, including the jackhmmer subsequence
        suffix.
  """
  name_to_sequence = collections.OrderedDict()
  for line in stockholm_string.splitlines():
    line = line.strip()
    if not line or line.startswith(('#', '//')):
      continue
    name, sequence = line.split()
    if name not in name_to_sequence:
      name_to_sequence[name] = ''
    name_to_sequence[name] += sequence

  msa = []
  deletion_matrix = []

  query = ''
  keep_columns = []
  for seq_index, sequence in enumerate(name_to_sequence.values()):
    if seq_index == 0:
      # Gather the columns with gaps from the query
      query = sequence
      keep_columns = [i for i, res in enumerate(query) if res != '-']

    # Remove the columns with gaps in the query from all sequences.
    aligned_sequence = ''.join([sequence[c] for c in keep_columns])

    msa.append(aligned_sequence)

    # Count the number of deletions w.r.t. query.
    deletion_vec = []
    deletion_count = 0
    for seq_res, query_res in zip(sequence, query):
      if seq_res != '-' or query_res != '-':
        if query_res == '-':
          deletion_count += 1
        else:
          deletion_vec.append(deletion_count)
          deletion_count = 0
    deletion_matrix.append(deletion_vec)

  return Msa(sequences=msa,
             deletion_matrix=deletion_matrix,
             descriptions=list(name_to_sequence.keys()))


def parse_a3m(a3m_string: str) -> Msa:
  """Parses sequences and deletion matrix from a3m format alignment.

  Args:
    a3m_string: The string contents of a a3m file. The first sequence in the
      file should be the query sequence.

  Returns:
    A tuple of:
      * A list of sequences that have been aligned to the query. These
        might contain duplicates.
      * The deletion matrix for the alignment as a list of lists. The element
        at `deletion_matrix[i][j]` is the number of residues deleted from
        the aligned sequence i at residue position j.
      * A list of descriptions, one per sequence, from the a3m file.
  """
  sequences, descriptions = parse_fasta(a3m_string)
  deletion_matrix = []
  for msa_sequence in sequences:
    deletion_vec = []
    deletion_count = 0
    for j in msa_sequence:
      if j.islower():
        deletion_count += 1
      else:
        deletion_vec.append(deletion_count)
        deletion_count = 0
    deletion_matrix.append(deletion_vec)

  # Make the MSA matrix out of aligned (deletion-free) sequences.
  deletion_table = str.maketrans('', '', string.ascii_lowercase)
  aligned_sequences = [s.translate(deletion_table) for s in sequences]
  return Msa(sequences=aligned_sequences,
             deletion_matrix=deletion_matrix,
             descriptions=descriptions)


def _convert_sto_seq_to_a3m(
    query_non_gaps: Sequence[bool], sto_seq: str) -> Iterable[str]:
  for is_query_res_non_gap, sequence_res in zip(query_non_gaps, sto_seq):
    if is_query_res_non_gap:
      yield sequence_res
    elif sequence_res != '-':
      yield sequence_res.lower()


def convert_stockholm_to_a3m(stockholm_format: str,
                             max_sequences: Optional[int] = None,
                             remove_first_row_gaps: bool = True) -> str:
  """Converts MSA in Stockholm format to the A3M format."""
  descriptions = {}
  sequences = {}
  reached_max_sequences = False

  for line in stockholm_format.splitlines():
    reached_max_sequences = max_sequences and len(sequences) >= max_sequences
    if line.strip() and not line.startswith(('#', '//')):
      # Ignore blank lines, markup and end symbols - remainder are alignment
      # sequence parts.
      seqname, aligned_seq = line.split(maxsplit=1)
      if seqname not in sequences:
        if reached_max_sequences:
          continue
        sequences[seqname] = ''
      sequences[seqname] += aligned_seq

  for line in stockholm_format.splitlines():
    if line[:4] == '#=GS':
      # Description row - example format is:
      # #=GS UniRef90_Q9H5Z4/4-78            DE [subseq from] cDNA: FLJ22755 ...
      columns = line.split(maxsplit=3)
      seqname, feature = columns[1:3]
      value = columns[3] if len(columns) == 4 else ''
      if feature != 'DE':
        continue
      if reached_max_sequences and seqname not in sequences:
        continue
      descriptions[seqname] = value
      if len(descriptions) == len(sequences):
        break

  # Convert sto format to a3m line by line
  a3m_sequences = {}
  if remove_first_row_gaps:
    # query_sequence is assumed to be the first sequence
    query_sequence = next(iter(sequences.values()))
    query_non_gaps = [res != '-' for res in query_sequence]
  for seqname, sto_sequence in sequences.items():
    # Dots are optional in a3m format and are commonly removed.
    out_sequence = sto_sequence.replace('.', '')
    if remove_first_row_gaps:
      out_sequence = ''.join(
          _convert_sto_seq_to_a3m(query_non_gaps, out_sequence))
    a3m_sequences[seqname] = out_sequence

  fasta_chunks = (f">{k} {descriptions.get(k, '')}\n{a3m_sequences[k]}"
                  for k in a3m_sequences)
  return '\n'.join(fasta_chunks) + '\n'  # Include terminating newline.


def _keep_line(line: str, seqnames: Set[str]) -> bool:
  """Function to decide which lines to keep."""
  if not line.strip():
    return True
  if line.strip() == '//':  # End tag
    return True
  if line.startswith('# STOCKHOLM'):  # Start tag
    return True
  if line.startswith('#=GC RF'):  # Reference Annotation Line
    return True
  if line[:4] == '#=GS':  # Description lines - keep if sequence in list.
    _, seqname, _ = line.split(maxsplit=2)
    return seqname in seqnames
  elif line.startswith('#'):  # Other markup - filter out
    return False
  else:  # Alignment data - keep if sequence in list.
    seqname = line.partition(' ')[0]
    return seqname in seqnames


def truncate_stockholm_msa(stockholm_msa_path: str, max_sequences: int) -> str:
  """Reads + truncates a Stockholm file while preventing excessive RAM usage."""
  seqnames = set()
  filtered_lines = []

  with open(stockholm_msa_path) as f:
    for line in f:
      if line.strip() and not line.startswith(('#', '//')):
        # Ignore blank lines, markup and end symbols - remainder are alignment
        # sequence parts.
        seqname = line.partition(' ')[0]
        seqnames.add(seqname)
        if len(seqnames) >= max_sequences:
          break

    f.seek(0)
    for line in f:
      if _keep_line(line, seqnames):
        filtered_lines.append(line)

  return ''.join(filtered_lines)


def remove_empty_columns_from_stockholm_msa(stockholm_msa: str) -> str:
  """Removes empty columns (dashes-only) from a Stockholm MSA."""
  processed_lines = {}
  unprocessed_lines = {}
  for i, line in enumerate(stockholm_msa.splitlines()):
    if line.startswith('#=GC RF'):
      reference_annotation_i = i
      reference_annotation_line = line
      # Reached the end of this chunk of the alignment. Process chunk.
      _, _, first_alignment = line.rpartition(' ')
      mask = []
      for j in range(len(first_alignment)):
        for _, unprocessed_line in unprocessed_lines.items():
          prefix, _, alignment = unprocessed_line.rpartition(' ')
          if alignment[j] != '-':
            mask.append(True)
            break
        else:  # Every row contained a hyphen - empty column.
          mask.append(False)
      # Add reference annotation for processing with mask.
      unprocessed_lines[reference_annotation_i] = reference_annotation_line

      if not any(mask):  # All columns were empty. Output empty lines for chunk.
        for line_index in unprocessed_lines:
          processed_lines[line_index] = ''
      else:
        for line_index, unprocessed_line in unprocessed_lines.items():
          prefix, _, alignment = unprocessed_line.rpartition(' ')
          masked_alignment = ''.join(itertools.compress(alignment, mask))
          processed_lines[line_index] = f'{prefix} {masked_alignment}'

      # Clear raw_alignments.
      unprocessed_lines = {}
    elif line.strip() and not line.startswith(('#', '//')):
      unprocessed_lines[i] = line
    else:
      processed_lines[i] = line
  return '\n'.join((processed_lines[i] for i in range(len(processed_lines))))


def deduplicate_stockholm_msa(stockholm_msa: str) -> str:
  """Remove duplicate sequences (ignoring insertions wrt query)."""
  sequence_dict = collections.defaultdict(str)

  # First we must extract all sequences from the MSA.
  for line in stockholm_msa.splitlines():
    # Only consider the alignments - ignore reference annotation, empty lines,
    # descriptions or markup.
    if line.strip() and not line.startswith(('#', '//')):
      line = line.strip()
      seqname, alignment = line.split()
      sequence_dict[seqname] += alignment

  seen_sequences = set()
  seqnames = set()
  # First alignment is the query.
  query_align = next(iter(sequence_dict.values()))
  mask = [c != '-' for c in query_align]  # Mask is False for insertions.
  for seqname, alignment in sequence_dict.items():
    # Apply mask to remove all insertions from the string.
    masked_alignment = ''.join(itertools.compress(alignment, mask))
    if masked_alignment in seen_sequences:
      continue
    else:
      seen_sequences.add(masked_alignment)
      seqnames.add(seqname)

  filtered_lines = []
  for line in stockholm_msa.splitlines():
    if _keep_line(line, seqnames):
      filtered_lines.append(line)

  return '\n'.join(filtered_lines) + '\n'


def _get_hhr_line_regex_groups(
    regex_pattern: str, line: str) -> Sequence[Optional[str]]:
  match = re.match(regex_pattern, line)
  if match is None:
    raise RuntimeError(f'Could not parse query line {line}')
  return match.groups()


def _update_hhr_residue_indices_list(
    sequence: str, start_index: int, indices_list: List[int]):
  """Computes the relative indices for each residue with respect to the original sequence."""
  counter = start_index
  for symbol in sequence:
    if symbol == '-':
      indices_list.append(-1)
    else:
      indices_list.append(counter)
      counter += 1


def _parse_hhr_hit(detailed_lines: Sequence[str]) -> TemplateHit:
  """Parses the detailed HMM HMM comparison section for a single Hit.

  This works on .hhr files generated from both HHBlits and HHSearch.

  Args:
    detailed_lines: A list of lines from a single comparison section between 2
      sequences (which each have their own HMM's)

  Returns:
    A dictionary with the information from that detailed comparison section

  Raises:
    RuntimeError: If a certain line cannot be processed
  """
  # Parse first 2 lines.
  number_of_hit = int(detailed_lines[0].split()[-1])
  name_hit = detailed_lines[1][1:]

  # Parse the summary line.
  pattern = (
      'Probab=(.*)[\t ]*E-value=(.*)[\t ]*Score=(.*)[\t ]*Aligned_cols=(.*)[\t'
      ' ]*Identities=(.*)%[\t ]*Similarity=(.*)[\t ]*Sum_probs=(.*)[\t '
      ']*Template_Neff=(.*)')
  match = re.match(pattern, detailed_lines[2])
  if match is None:
    raise RuntimeError(
        'Could not parse section: %s. Expected this: \n%s to contain summary.' %
        (detailed_lines, detailed_lines[2]))
  (_, _, _, aligned_cols, _, _, sum_probs, _) = [float(x)
                                                 for x in match.groups()]

  # The next section reads the detailed comparisons. These are in a 'human
  # readable' format which has a fixed length. The strategy employed is to
  # assume that each block starts with the query sequence line, and to parse
  # that with a regexp in order to deduce the fixed length used for that block.
  query = ''
  hit_sequence = ''
  indices_query = []
  indices_hit = []
  length_block = None

  for line in detailed_lines[3:]:
    # Parse the query sequence line
    if (line.startswith('Q ') and not line.startswith('Q ss_dssp') and
        not line.startswith('Q ss_pred') and
        not line.startswith('Q Consensus')):
      # Thus the first 17 characters must be 'Q <query_name> ', and we can parse
      # everything after that.
      #              start    sequence       end       total_sequence_length
      patt = r'[\t ]*([0-9]*) ([A-Z-]*)[\t ]*([0-9]*) \([0-9]*\)'
      groups = _get_hhr_line_regex_groups(patt, line[17:])

      # Get the length of the parsed block using the start and finish indices,
      # and ensure it is the same as the actual block length.
      start = int(groups[0]) - 1  # Make index zero based.
      delta_query = groups[1]
      end = int(groups[2])
      num_insertions = len([x for x in delta_query if x == '-'])
      length_block = end - start + num_insertions
      assert length_block == len(delta_query)

      # Update the query sequence and indices list.
      query += delta_query
      _update_hhr_residue_indices_list(delta_query, start, indices_query)

    elif line.startswith('T '):
      # Parse the hit sequence.
      if (not line.startswith('T ss_dssp') and
          not line.startswith('T ss_pred') and
          not line.startswith('T Consensus')):
        # Thus the first 17 characters must be 'T <hit_name> ', and we can
        # parse everything after that.
        #              start    sequence       end     total_sequence_length
        patt = r'[\t ]*([0-9]*) ([A-Z-]*)[\t ]*[0-9]* \([0-9]*\)'
        groups = _get_hhr_line_regex_groups(patt, line[17:])
        start = int(groups[0]) - 1  # Make index zero based.
        delta_hit_sequence = groups[1]
        assert length_block == len(delta_hit_sequence)

        # Update the hit sequence and indices list.
        hit_sequence += delta_hit_sequence
        _update_hhr_residue_indices_list(delta_hit_sequence, start, indices_hit)

  return TemplateHit(
      index=number_of_hit,
      name=name_hit,
      aligned_cols=int(aligned_cols),
      sum_probs=sum_probs,
      query=query,
      hit_sequence=hit_sequence,
      indices_query=indices_query,
      indices_hit=indices_hit,
  )


def parse_hhr(hhr_string: str) -> Sequence[TemplateHit]:
  """Parses the content of an entire HHR file."""
  lines = hhr_string.splitlines()

  # Each .hhr file starts with a results table, then has a sequence of hit
  # "paragraphs", each paragraph starting with a line 'No <hit number>'. We
  # iterate through each paragraph to parse each hit.

  block_starts = [i for i, line in enumerate(lines) if line.startswith('No ')]

  hits = []
  if block_starts:
    block_starts.append(len(lines))  # Add the end of the final block.
    for i in range(len(block_starts) - 1):
      hits.append(_parse_hhr_hit(lines[block_starts[i]:block_starts[i + 1]]))
  return hits


def parse_e_values_from_tblout(tblout: str) -> Dict[str, float]:
  """Parse target to e-value mapping parsed from Jackhmmer tblout string."""
  e_values = {'query': 0}
  lines = [line for line in tblout.splitlines() if line[0] != '#']
  # As per http://eddylab.org/software/hmmer/Userguide.pdf fields are
  # space-delimited. Relevant fields are (1) target name:  and
  # (5) E-value (full sequence) (numbering from 1).
  for line in lines:
    fields = line.split()
    e_value = fields[4]
    target_name = fields[0]
    e_values[target_name] = float(e_value)
  return e_values


def _get_indices(sequence: str, start: int) -> List[int]:
  """Returns indices for non-gap/insert residues starting at the given index."""
  indices = []
  counter = start
  for symbol in sequence:
    # Skip gaps but add a placeholder so that the alignment is preserved.
    if symbol == '-':
      indices.append(-1)
    # Skip deleted residues, but increase the counter.
    elif symbol.islower():
      counter += 1
    # Normal aligned residue. Increase the counter and append to indices.
    else:
      indices.append(counter)
      counter += 1
  return indices


@dataclasses.dataclass(frozen=True)
class HitMetadata:
  pdb_id: str
  chain: str
  start: int
  end: int
  length: int
  text: str


def _parse_hmmsearch_description(description: str) -> HitMetadata:
  """Parses the hmmsearch A3M sequence description line."""
  # Example 1: >4pqx_A/2-217 [subseq from] mol:protein length:217  Free text
  # Example 2: >5g3r_A/1-55 [subseq from] mol:protein length:352
  match = re.match(
      r'^>?([a-z0-9]+)_(\w+)/([0-9]+)-([0-9]+).*protein length:([0-9]+) *(.*)$',
      description.strip())

  if not match:
    raise ValueError(f'Could not parse description: "{description}".')

  return HitMetadata(
      pdb_id=match[1],
      chain=match[2],
      start=int(match[3]),
      end=int(match[4]),
      length=int(match[5]),
      text=match[6])


def parse_hmmsearch_a3m(query_sequence: str,
                        a3m_string: str,
                        skip_first: bool = True) -> Sequence[TemplateHit]:
  """Parses an a3m string produced by hmmsearch.

  Args:
    query_sequence: The query sequence.
    a3m_string: The a3m string produced by hmmsearch.
    skip_first: Whether to skip the first sequence in the a3m string.

  Returns:
    A sequence of `TemplateHit` results.
  """
  # Zip the descriptions and MSAs together, skip the first query sequence.
  parsed_a3m = list(zip(*parse_fasta(a3m_string)))
  if skip_first:
    parsed_a3m = parsed_a3m[1:]

  indices_query = _get_indices(query_sequence, start=0)

  hits = []
  for i, (hit_sequence, hit_description) in enumerate(parsed_a3m, start=1):
    if 'mol:protein' not in hit_description:
      continue  # Skip non-protein chains.
    metadata = _parse_hmmsearch_description(hit_description)
    # Aligned columns are only the match states.
    aligned_cols = sum([r.isupper() and r != '-' for r in hit_sequence])
    indices_hit = _get_indices(hit_sequence, start=metadata.start - 1)

    hit = TemplateHit(
        index=i,
        name=f'{metadata.pdb_id}_{metadata.chain}',
        aligned_cols=aligned_cols,
        sum_probs=None,
        query=query_sequence,
        hit_sequence=hit_sequence.upper(),
        indices_query=indices_query,
        indices_hit=indices_hit,
    )
    hits.append(hit)

  return hits
