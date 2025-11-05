import subprocess
import os
from bs4 import BeautifulSoup

def create_config_file(srr_id, sequence_length, config_path='config.txt'):
    config_content = f"""Project:
-----------------------
Project name          = {srr_id}_mito
Type                  = mito
Genome Range          = 30000-55000
K-mer                 = 25
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = /home/labpc1/seed.fasta
Extend seed directly  = no
Reference sequence    = /home/labpc1/L757_Mito_Ref.fasta
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = {sequence_length}
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = /home/labpc1/{srr_id}_1.fastq
Reverse reads         = /home/labpc1/{srr_id}_2.fastq
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
Output path           = /media/labpc1/MARIANA_2/Q30/Candida_Assembly/{srr_id}/NOVOPlasty/
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

        # Check if the length contains a hyphen
        if '-' in sequence_length:
            # Split the string and take the second part
            sequence_length = sequence_length.split('-')[1].strip()

        sequence_length = int(sequence_length)  # Convert to integer for FLASH
    else:
        sequence_length = None

    return sequence_length

def log_skipped_srr(srr_id, reason, log_file='skipped_srr_ids.txt'):
    with open(log_file, 'a') as file:
        file.write(f"{srr_id} skipped: {reason}\n")

def main():
    # Read SRR IDs from the text file
    with open('SRR_list.txt', 'r') as file:
        srr_ids = [line.strip() for line in file if line.strip()]  # Read non-empty lines

    for srr_id in srr_ids:
        print(f"Processing SRR ID: {srr_id}")

        # Run prefetch command
        result = subprocess.run(['prefetch', srr_id])
        if result.returncode != 0:
            print(f"Error downloading {srr_id} with prefetch. Skipping to next SRR ID.")
            log_skipped_srr(srr_id, "prefetch failed")
            continue

        # Run fastq-dump command
        result = subprocess.run(['fastq-dump', '--split-3', srr_id])
        if result.returncode != 0 or not (os.path.exists(f'{srr_id}_1.fastq') and os.path.exists(f'{srr_id}_2.fastq')):
            print(f"Error running fastq-dump for {srr_id}. Skipping to next SRR ID.")
            log_skipped_srr(srr_id, "fastq-dump failed")
            continue

        # Run FastQC command
        result = subprocess.run(['fastqc', '-t', '10', f'{srr_id}_1.fastq', f'{srr_id}_2.fastq'])
        if result.returncode != 0 or not os.path.exists(f'{srr_id}_1_fastqc.html'):
            print(f"Error running FastQC for {srr_id}. Skipping to next SRR ID.")
            log_skipped_srr(srr_id, "FastQC failed")
            continue

        # Extract sequence length from the FastQC HTML file
        sequence_length = extract_sequence_length(f'{srr_id}_1_fastqc.html')
        if sequence_length is None:
            print(f"Failed to extract sequence length for {srr_id}. Skipping to next SRR ID.")
            log_skipped_srr(srr_id, "Failed to extract sequence length")
            continue

        # Run TrimGalore
        result = subprocess.run(['trim_galore', '--paired', '-q', '30', '--length', '20', '-j', '4', f'{srr_id}_1.fastq', f'{srr_id}_2.fastq'])
        if result.returncode != 0:
            print(f"Error running TrimGalore for {srr_id}. Skipping to next SRR ID.")
            log_skipped_srr(srr_id, "TrimGalore failed")
            continue

        # Run FLASH
        result = subprocess.run(['flash', '-m', '25', '-M', str(sequence_length), '-o', f'flash_{srr_id}', f'{srr_id}_1_val_1.fq', f'{srr_id}_2_val_2.fq'])
        if result.returncode != 0:
            print(f"Error running FLASH for {srr_id}. Skipping to next SRR ID.")
            log_skipped_srr(srr_id, "FLASH failed")
            continue

        # Run SPAdes
        result = subprocess.run(['spades.py', '-1', f'flash_{srr_id}.notCombined_1.fastq', '-2', f'flash_{srr_id}.notCombined_2.fastq', '--merged', f'flash_{srr_id}.extendedFrags.fastq', '-o', f'{srr_id}_q30_genome'])
        if result.returncode != 0:
            print(f"Error running SPAdes for {srr_id}. Skipping to next SRR ID.")
            log_skipped_srr(srr_id, "SPAdes failed")
            continue

        # Run NOVOPlasty
        # Create config file for NOVOPlasty
        # Create directories if they don't exist
        output_directory_1 = f'/media/labpc1/MARIANA_2/Q30/Candida_Assembly/{srr_id}/NOVOPlasty/'
        os.makedirs(output_directory_1, exist_ok=True)
        create_config_file(srr_id, sequence_length)

        result = subprocess.run(['NOVOPlasty4.3.5.pl', '-c', 'config.txt'])
        if result.returncode != 0:
            print(f"Error running NOVOPlasty for {srr_id}. Skipping to next SRR ID.")
            log_skipped_srr(srr_id, "NOVOPlasty failed")
            continue

        # Remove unnecessary files
        files_to_remove = [
            f'{srr_id}_1.fastq', f'{srr_id}_2.fastq', f'flash_{srr_id}.notCombined_1.fastq',
            f'flash_{srr_id}.notCombined_2.fastq', f'{srr_id}_1_val_1.fq', f'{srr_id}_2_val_2.fq',
            f'flash_{srr_id}.extendedFrags.fastq', f'flash_{srr_id}.hist', f'flash_{srr_id}.histogram',
            f'sra-toolkit/sra/{srr_id}.sra', f'{srr_id}_1_fastqc.zip', f'{srr_id}_2_fastqc.zip'
        ]
        for file in files_to_remove:
            subprocess.run(['rm', file])

        # Create directories if they don't exist
        output_directory_2 = f'/media/labpc1/MARIANA_2/Q30/Candida_Assembly/{srr_id}/'
        os.makedirs(output_directory_2, exist_ok=True)

        # Move files to MARIANA_2 external drive
        files_to_move = [
            f'{srr_id}_q30_genome', f'{srr_id}_1.fastq_trimming_report.txt',
            f'{srr_id}_2.fastq_trimming_report.txt', f'{srr_id}_1_fastqc.html',
            f'{srr_id}_2_fastqc.html', f'{srr_id}_GenoSplit'
        ]
        for file in files_to_move:
            subprocess.run(['mv', file, output_directory_2])

if __name__ == "__main__":
    main()

