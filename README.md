# Tiger Phylogenetics
**Tugas Besar IF3211 Komputasi Domain Spesifik**

## Anggota Kelompok

| NIM | Nama |
|---|---|
| 13523007 | Ranashahira Reztaputri |
| 13523043 | Najwa Kahani Fatima |
| 13523050 | Mayla Yaffa Ludmilla |
| 13523058 | Noumisyifa Nabila Nareswari |

Membangun pohon filogenetik harimau Nusantara dan Asia menggunakan sekuens mitokondrial cytochrome b, dengan tujuan mengukur jarak genetik antar subspesies dan menganalisis gap evolusioner akibat kepunahan Harimau Jawa.

---

## Dataset

Sekuens diambil dari NCBI GenBank menggunakan accession ID berikut:

| Label | Subspesies | Accession | Status |
|---|---|---|---|
| P_t_sumatrae | Harimau Sumatera | AF053054.1 | Kritis (CR) |
| P_t_sondaica | Harimau Jawa | OQ601561.1 | **Punah** |
| P_t_tigris | Harimau Bengal | AF053053.1 | Endangered |
| P_t_altaica | Harimau Amur/Siberia | AF053051.1 | Endangered |
| P_t_corbetti | Harimau Indochinese | AF053050.1 | Endangered |
| P_t_jacksoni | Harimau Malayan | EU184702.1 | Kritis (CR) |

---

## Pipeline

```
fetch_sequence.py   →  smith_waterman.py  →  run_mafft.py
(Fetch dari NCBI)      (Pairwise alignment)   (Multiple alignment)
                                                      ↓
visualizer.py       ←  evolutionary_gap.py ←  kimura_distance.py
(Pohon + plot)         (Gap analysis)         (K2P distance matrix)
                              ↑
                       bootstrap_sampling.py  ←  neighbor_joining.py
                       (Validasi 1000x)           (Bangun pohon NJ)
```

| # | File | Fungsi |
|---|---|---|
| 1 | `fetch_sequence.py` | Ambil sekuens FASTA dari NCBI via Biopython Entrez |
| 2 | `smith_waterman.py` | Pairwise alignment dengan algoritma Smith-Waterman |
| 3 | `run_mafft.py` | Multiple sequence alignment dengan MAFFT + trimming kolom gap |
| 4 | `kimura_distance.py` | Hitung jarak genetik Kimura 2-Parameter (K2P) |
| 5 | `neighbor_joining.py` | Bangun pohon filogenetik Neighbor-Joining dari scratch |
| 6 | `bootstrap_sampling.py` | Bootstrap 1000x untuk validasi confidence setiap cabang |
| 7 | `evolutionary_gap.py` | Analisis gap evolusioner dan distinctiveness tiap subspesies |
| 8 | `visualizer.py` | Render pohon filogenetik ke PNG dengan Matplotlib |

---

## Cara Menjalankan

### 1. Clone repo
```bash
git clone https://github.com/<username>/tiger-phylogenetics.git
cd tiger-phylogenetics
```

### 2. Install dependencies
```bash
pip install biopython matplotlib numpy scipy
```

### 3. Install MAFFT (Windows)
- Download [mafft-7.526-win64-signed](https://mafft.cbrc.jp/alignment/software/windows.html)
- Ekstrak ke direktori pilihan kamu
- Buka `run_mafft.py` dan sesuaikan `MAFFT_PATH` dengan path hasil ekstrak:
  ```python
  MAFFT_PATH = r"C:\path\ke\mafft-win\mafft.bat"
  ```

> **Linux/Mac:** Install via conda: `conda install -c bioconda mafft`

### 4. Jalankan program
```bash
python main.py
```

Output akan tersimpan di:
- `data/` — sekuens mentah, hasil alignment, distance matrix
- `results/phylogenetic_tree.png` — visualisasi pohon filogenetik

---

## Struktur Folder

```
tiger-phylogenetics/
├── data/
│   ├── accession_ids.json          # Mapping label -> accession NCBI
│   ├── raw_sequences.fasta         # Sekuens mentah hasil fetch
│   ├── aligned_raw.fasta           # Output MAFFT sebelum trimming
│   ├── aligned_sequences.fasta     # Alignment final (setelah trim gap >50%)
│   ├── pairwise_scores.json        # Skor Smith-Waterman tiap pasang
│   └── pairwise_alignments.txt     # Detail alignment tiap pasang
├── results/
│   └── phylogenetic_tree.png       # Output visualisasi pohon
├── fetch_sequence.py
├── smith_waterman.py
├── run_mafft.py
├── kimura_distance.py
├── neighbor_joining.py
├── bootstrap_sampling.py
├── evolutionary_gap.py
├── visualizer.py
├── main.py
└── README.md
```

---

## Output

**Pohon filogenetik** dengan kode warna:
- 🟢 Hijau — Harimau Sumatera
- 🔴 Merah — Harimau Jawa (punah)
- ⚫ Abu-abu — Subspesies Asia lainnya

Angka di tiap cabang adalah **bootstrap support** dari 1000 ulangan:
- **Biru ≥70%** — dukungan kuat
- **Kuning <70%** — dukungan lemah

**Evolutionary gap analysis** mencakup:
- Branch length yang hilang akibat kepunahan Harimau Jawa
- Proporsi total evolutionary loss dari keseluruhan pohon
- Evolutionary distinctiveness score tiap subspesies

---

## Dependencies

| Library | Versi minimum | Kegunaan |
|---|---|---|
| biopython | ≥1.79 | Fetch NCBI, parsing FASTA, Phylo |
| numpy | ≥1.21 | Matriks jarak, komputasi K2P |
| matplotlib | ≥3.4 | Visualisasi pohon |
| MAFFT | 7.x | Multiple sequence alignment (eksternal) |

---

## Referensi

- Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. *Journal of Molecular Evolution*, 16, 111–120.
- Saitou, N., & Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution*, 4(4), 406–425.
- Smith, T.F., & Waterman, M.S. (1981). Identification of common molecular subsequences. *Journal of Molecular Biology*, 147(1), 195–197.
- NCBI GenBank: https://www.ncbi.nlm.nih.gov/nucleotide/
