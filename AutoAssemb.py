import argparse
import subprocess
import os
from bs4 import BeautifulSoup

def create_config_file(srr_id, sequence_length, config_path, seed_path, refseq_path, output_path):
    config_content = f"""Project:
-----------------------
Project name          = {srr_id}_mito
Type                  = mito
Genome Range          = 30000-55000
K-mer                 = 25
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = {seed_path}
Extend seed directly  = no
Reference sequence    = {refseq_path}
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = {sequence_length}
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = {os.path.join(output_path, srr_id)}_1.fastq
Reverse reads         = {os.path.join(output_path, srr_id)}_2.fastq
Store Hash            =

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
Reduce ambigious N's  =
Output path           = {os.path.join(output_path, srr_id, 'NOVOPlasty')}
"""
    with open(config_path, 'w') as file:
        file.write(config_content)

def extract_sequence_length(html_file):
    with open(html_file, 'r') as file:
        html_content = file.read()

    soup = BeautifulSoup(html_content, 'html.parser')
    sequence_length_tag = soup.find('td', string='Sequence length')

    if sequence_length_tag:
        length_tag = sequence_length_tag.find_next('td')
        sequence_length = length_tag.text.strip()

        # Handle hyphen cases
        if '-' in sequence_length:
            sequence_length = sequence_length.split('-')[1].strip()

        sequence_length = int(sequence_length)
    else:
        sequence_length = None

    return sequence_length

def log_skipped_srr(srr_id, reason, log_file):
    with open(log_file, 'a') as file:
        file.write(f"{srr_id} skipped: {reason}\n")

def main():
    parser = argparse.ArgumentParser(description='AutoAssemb user-friendly pipeline')
    parser.add_argument('--srr-list', required=True, help='File with SRR IDs and isolate names (one per line, comma separated: SRR,ID)')
    parser.add_argument('--output-dir', required=True, help='Directory for output and final results')
    parser.add_argument('--work-dir', required=True, help='Working directory for temp/intermediate files')
    parser.add_argument('--seed-file', required=True, help='Path to seed FASTA')
    parser.add_argument('--refseq-file', required=True, help='Path to reference FASTA')
    parser.add_argument('--log-file', default='skipped_srr_ids.txt', help='Log file for skipped SRRs')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.work_dir, exist_ok=True)

    with open(args.srr_list, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    srr_entries = [line.split(',') for line in lines]

    for srr_id, isolate_id in srr_entries:
        print(f"\nProcessing: SRR ID: {srr_id}, Isolate ID: {isolate_id}")

        # All temp files go to the working dir!
        wdir = args.work_dir
        outdir = os.path.join(args.output_dir, isolate_id)
        novoplasty_dir = os.path.join(outdir, "NOVOPlasty")
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(novoplasty_dir, exist_ok=True)

        fq1 = os.path.join(wdir, f'{srr_id}_1.fastq')
        fq2 = os.path.join(wdir, f'{srr_id}_2.fastq')

        # 1. Prefetch, 2. fastq-dump, 3. FastQC
        result = subprocess.run(['prefetch', srr_id], cwd=wdir)
        if result.returncode != 0:
            print(f"Error prefetch {srr_id}")
            log_skipped_srr(srr_id, f"prefetch failed", args.log_file)
            continue

        result = subprocess.run(['fastq-dump', '--split-3', srr_id], cwd=wdir)
        if result.returncode != 0 or not (os.path.exists(fq1) and os.path.exists(fq2)):
            print(f"Error fastq-dump {srr_id}")
            log_skipped_srr(srr_id, f"fastq-dump failed", args.log_file)
            continue

        result = subprocess.run(['fastqc', '-t', '10', fq1, fq2], cwd=wdir)
        fq1_fastqc = os.path.join(wdir, f'{srr_id}_1_fastqc.html')
        if result.returncode != 0 or not os.path.exists(fq1_fastqc):
            print(f"Error FastQC {srr_id}")
            log_skipped_srr(srr_id, f"FastQC failed", args.log_file)
            continue

        sequence_length = extract_sequence_length(fq1_fastqc)
        if sequence_length is None:
            print(f"Failed to extract sequence length {srr_id}")
            log_skipped_srr(srr_id, f"Failed to extract sequence length", args.log_file)
            continue

        # Run TrimGalore
        result = subprocess.run(['trim_galore', '--paired', '-q', '30', '--length', '20', '-j', '4', fq1, fq2], cwd=wdir)
        fq1_val = os.path.join(wdir, f'{srr_id}_1_val_1.fq')
        fq2_val = os.path.join(wdir, f'{srr_id}_2_val_2.fq')
        if result.returncode != 0:
            print(f"Error TrimGalore {srr_id}")
            log_skipped_srr(srr_id, f"TrimGalore failed", args.log_file)
            continue

        # FLASH
        result = subprocess.run(['flash', '-m', '25', '-M', str(sequence_length), '-o', f'flash_{srr_id}', fq1_val, fq2_val], cwd=wdir)
        if result.returncode != 0:
            print(f"Error FLASH {srr_id}")
            log_skipped_srr(srr_id, f"FLASH failed", args.log_file)
            continue

        # SPAdes
        result = subprocess.run(
            ['spades.py', '-1', f'flash_{srr_id}.notCombined_1.fastq', '-2', f'flash_{srr_id}.notCombined_2.fastq', '--merged', f'flash_{srr_id}.extendedFrags.fastq', '-o', f'spades_{srr_id}'],
            cwd=wdir
        )
        if result.returncode != 0:
            print(f"Error SPAdes {srr_id}")
            log_skipped_srr(srr_id, f"SPAdes failed", args.log_file)
            continue

        # NOVOPlasty config and run
        config_path = os.path.join(wdir, f'config_{srr_id}.txt')
        create_config_file(srr_id, sequence_length, config_path, args.seed_file, args.refseq_file, outdir)
        result = subprocess.run(['NOVOPlasty4.3.5.pl', '-c', config_path], cwd=wdir)
        if result.returncode != 0:
            print(f"Error NOVOPlasty {srr_id}")
            log_skipped_srr(srr_id, f"NOVOPlasty failed", args.log_file)
            continue

        # Move outputs to output directory
        files_to_copy = [
            f'{srr_id}_q30_genome',
            f'{srr_id}_1.fastq_trimming_report.txt',
            f'{srr_id}_2.fastq_trimming_report.txt',
            f'{srr_id}_1_fastqc.html',
            f'{srr_id}_2_fastqc.html',
            # Add more outputs as needed
        ]
        for file in files_to_copy:
            src = os.path.join(wdir, file)
            dst = os.path.join(outdir, file)
            if os.path.exists(src):
                os.rename(src, dst)
        print(f"Finished {srr_id}")

if __name__ == '__main__':
    main()
