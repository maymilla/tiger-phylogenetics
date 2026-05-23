import numpy as np
import json
from Bio import SeqIO
from itertools import combinations

MATCH =  2
MISMATCH = -1
GAP = -2

def smith_waterman(seq1, seq2):
  n = len(seq1) + 1
  m = len(seq2) + 1
  
  matrix = np.zeros((n, m), dtype=int)
  
  for i in range(1, n):
    for j in range(1, m):
      diagonal = matrix[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH)
      up = matrix[i-1][j] + GAP 
      left = matrix[i][j-1] + GAP 
      
      matrix[i][j] = max(0, diagonal, up, left)
  
  return matrix
  
def traceback(matrix, seq1, seq2):
  score = np.max(matrix)
  pos = np.unravel_index(np.argmax(matrix), matrix.shape)
  i, j = pos
  
  aligned1 = ""
  aligned2 = ""
  
  while matrix[i][j] != 0:
    if seq1[i-1] == seq2[j-1]:
      aligned1 = seq1[i-1] + aligned1
      aligned2 = seq2[j-1] + aligned2
      i -= 1
      j -= 1
    elif matrix[i-1][j] + GAP == matrix[i][j]:
      aligned1 = seq1[i-1] + aligned1
      aligned2 = "-" + aligned2
      i -= 1
    else:
      aligned1 = "-" + aligned1
      aligned2 = seq2[j-1] + aligned2
      j -= 1
  
  return aligned1, aligned2, int(score), pos

def main():
  fasta_path = "data/raw_sequences.fasta"
  records = list(SeqIO.parse(fasta_path, "fasta"))
  
  print(f"Loaded {len(records)} sekuens")
  
  pairwise_scores = {}   
  alignment_lines = []   
  
  pairs = list(combinations(records, 2))
  total = len(pairs)
  
  for idx, (rec1, rec2) in enumerate(pairs, 1):
    label1 = rec1.id
    label2 = rec2.id
    seq1 = str(rec1.seq)
    seq2 = str(rec2.seq)
    
    print(f"[{idx}/{total}] Aligning {label1} vs {label2} ...")
    
    matrix_sw = smith_waterman(seq1, seq2)
    aligned1, aligned2, score, end_pos = traceback(matrix_sw, seq1, seq2)
    
    matches = sum(a == b for a, b in zip(aligned1, aligned2) if a != '-' and b != '-')
    aln_len = len(aligned1)
    identity = round(matches / aln_len * 100, 2) if aln_len > 0 else 0
    
    pair_key = f"{label1}_vs_{label2}"
    pairwise_scores[pair_key] = {
      "score":        score,
      "identity_pct": identity,
      "aln_length":   aln_len,
      "end_pos_seq1": int(end_pos[0]),
      "end_pos_seq2": int(end_pos[1])
    }
    
    alignment_lines.append(f"### {pair_key}")
    alignment_lines.append(f"Score: {score} | Identity: {identity}% | Length: {aln_len}")
    alignment_lines.append(f"Seq1: {aligned1[:80]}...")  
    alignment_lines.append(f"Seq2: {aligned2[:80]}...")
    alignment_lines.append("")
    
  json_path = "data/pairwise_scores.json"
  with open(json_path, "w") as f:
    json.dump(pairwise_scores, f, indent=2)
  print(f"\nScores disimpan ke {json_path}")
  
  txt_path = "data/pairwise_alignments.txt"
  with open(txt_path, "w") as f:
    f.write("\n".join(alignment_lines))
  print(f"Alignments disimpan ke {txt_path}")

if __name__ == "__main__":
  main()