# AutoAssemb

**AutoAssemb** is a user-friendly, automated genome assembly pipeline designed for high-throughput sequencing projects. Given a simple input file with SRA Run accessions and isolate IDs, along with paths for your seed and reference FASTA files, AutoAssemb streamlines the process from raw data download to organelle genome assembly and report collation.

## Features

- **Automated download** of raw SRA data
- **Quality control, trimming, merging, and assembly** (using FastQC, TrimGalore, FLASH, and SPAdes)
- **Organelle assembly using NOVOPlasty**
- **Customizable paths** for input, output, seed, reference, and working directories
- **Error logging** for skipped/failed samples

## Requirements

- Python 3.x
- [sra-tools](https://github.com/ncbi/sra-tools) (`prefetch`, `fastq-dump`)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [FLASH](https://ccb.jhu.edu/software/FLASH/)
- [SPAdes](https://github.com/ablab/spades)
- [NOVOPlasty](https://github.com/ndierckx/NOVOPlasty)
- [BeautifulSoup4](https://pypi.org/project/beautifulsoup4/)

> **All tools must be installed and available in your system `$PATH`.**

Install Python dependencies using pip (if needed):

```sh
pip install beautifulsoup4
```

## Input Preparation

- Create a plain text input file (e.g., `srr_list.txt`) with each line as:  
  ```
  SRR123456,IsolateName
  SRR987654,AnotherIsolate
  ```

- Obtain your **seed FASTA** and **reference sequence FASTA**.

- Decide where you want the **output** (final results) and **working** (temporary/intermediate) directories.

## Usage

```sh
python AutoAssemb.py \
  --srr-list path/to/srr_list.txt \
  --output-dir path/to/final/output \
  --work-dir path/to/temp/working \
  --seed-file path/to/seed.fasta \
  --refseq-file path/to/reference.fasta
```
**Required arguments:**
- `--srr-list`      Path to your SRA accession and isolate list file
- `--output-dir`    Directory for all final output files
- `--work-dir`      Directory for temporary/intermediate files
- `--seed-file`     Path to your seed FASTA for NOVOPlasty
- `--refseq-file`   Reference sequence FASTA for NOVOPlasty

Optional:
- `--log-file`      Path for logging skipped samples (default: `skipped_srr_ids.txt`)

## Output

- Results and reports for each isolate will be saved in the structure:
  ```
  output-dir/
    IsolateName/
      NOVOPlasty/
      [other results and reports]
  ```
- Skipped/failed samples are logged for easy troubleshooting.

## Example

```sh
python AutoAssemb.py \
  --srr-list my_samples.csv \
  --output-dir /home/user/final_results \
  --work-dir /home/user/work_tmp \
  --seed-file /home/user/seed.fa \
  --refseq-file /home/user/refseq.fa
```

## Notes

- Project is robust: if one sample fails, processing continues without losing previous results.
- Compatible with any number of samples—just list them in your input file.

## License

[MIT](LICENSE)

----

*For any questions or suggestions, please open an issue or pull request!*
