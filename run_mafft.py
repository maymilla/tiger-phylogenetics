import subprocess
import sys
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

INPUT_FASTA = "data/raw_sequences.fasta"
RAW_ALIGNED = "data/aligned_raw.fasta"       
OUTPUT_FASTA = "data/aligned_sequences.fasta" 
MAFFT_PATH = r"C:\Users\raiha\Downloads\mafft-7.526-win64-signed\mafft-win\mafft.bat"
 

def run_mafft(input_path, output_path):
  print("Menjalankan MAFFT ...")
  
  result = subprocess.run(
    [MAFFT_PATH, "--auto", "--quiet", input_path],
    capture_output=True,
    text=True,
    shell=True  
  )
  
  if result.returncode != 0:
    print("[ERROR] MAFFT gagal:")
    print(result.stderr)
    sys.exit(1)
  
  with open(output_path, "w") as f:
    f.write(result.stdout)

  print(f"MAFFT selesai. Output: {output_path}")
  
# hapus kolom alignment yang >50% isinya gap. 
def trim_gappy_columns(input_path, output_path, gap_threshold=0.5):
  print(f"Trimming kolom dengan gap > {gap_threshold*100:.0f}% ...")
  
  alignment = list(SeqIO.parse(input_path, "fasta"))
  n_seq = len(alignment)
  aln_len = len(alignment[0].seq)
  
  print(f"Sebelum trim: {n_seq} sekuens, panjang {aln_len} bp")
  
  keep_cols = []
  for col_idx in range(aln_len):
    gap_count = sum(1 for rec in alignment if rec.seq[col_idx] == '-')
    gap_fraction = gap_count / n_seq
    
    if gap_fraction <= gap_threshold:
      keep_cols.append(col_idx)
  
  print(f"Kolom dipertahankan: {len(keep_cols)}/{aln_len}")
  
  trimmed_records = []
  for rec in alignment:
    trimmed_seq = "".join(rec.seq[i] for i in keep_cols)
    trimmed_records.append(SeqRecord(
      Seq(trimmed_seq),
      id=rec.id,
      description=rec.description
    ))
  
  SeqIO.write(trimmed_records, output_path, "fasta")
  print(f"Hasil trim disimpan. Output: {output_path}")
  
  return trimmed_records

# pastikan semua sekuens panjangnya sama setelah trim
def validate_alignment(records):
  print("[INFO] Validasi alignment ...")
  
  lengths = set(len(rec.seq) for rec in records)
  
  if len(lengths) != 1:
    print(f"[ERROR] Panjang sekuens tidak konsisten: {lengths}")
    sys.exit(1)
  
  final_len = lengths.pop()
  print(f"Semua {len(records)} sekuens panjang {final_len} bp — siap handoff")
  
  print("\nRingkasan:")
  for rec in records:
    gap_count = rec.seq.count('-')
    gap_pct = gap_count / len(rec.seq) * 100
    print(f"  {rec.id:<20} {len(rec.seq)} bp  |  gap: {gap_count} ({gap_pct:.1f}%)")
    
def main(input_path=INPUT_FASTA, raw_aligned_path=RAW_ALIGNED, output_path=OUTPUT_FASTA):
  if not os.path.exists(input_path):
    print(f"File tidak ditemukan: {input_path}")
    sys.exit(1)
  
  run_mafft(input_path, raw_aligned_path)
  
  trimmed = trim_gappy_columns(raw_aligned_path, output_path)
  
  validate_alignment(trimmed)
  
  print(f"\nHasil: {output_path}")
  return output_path

if __name__ == "__main__":
  main()