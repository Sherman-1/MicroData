import os
import subprocess
import sys
from pathlib import Path
import shutil
from Bio import SeqIO
import argparse


class MMseqs2:
    def __init__(self, threads, mmseqs2_path, cleanup=True):

        self.mmseqs2_path = mmseqs2_path
        self.threads = threads
        self.cleanup = cleanup
        self.dir = Path(os.getcwd()) / "mmseqs_tmp"
        self.dir.mkdir(parents=True, exist_ok=True)

        print(f"MMseqs running in {os.getcwd()}")

    def _run_command(self, command, log_file):

        print(f"Running command: {command}")
        with open(log_file, 'w') as log:
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            for stdout_line in iter(process.stdout.readline, ""):
                log.write(stdout_line)
            for stderr_line in iter(process.stderr.readline, ""):
                print(stderr_line, end="", file=sys.stderr)
                log.write(stderr_line)
            process.stdout.close()
            process.stderr.close()
            return_code = process.wait()
            if return_code != 0:
                raise Exception(f"Command '{command}' failed with return code {return_code}")

    def createdb(self, fasta, log_file='createDB.log'):
        path = Path(fasta)
        if not path.exists():
            raise FileNotFoundError(f"FASTA file {fasta} does not exist.")
        command = f"{self.mmseqs2_path} createdb {path} {self.dir}/db"
        self._run_command(command, log_file)

    def cluster(self, coverage, identity, cov_mode, log_file='clustering.log'):
        command = (
            f"{self.mmseqs2_path} cluster {self.dir}/db {self.dir}/db_clu {self.dir}/tmp -c {coverage} --min-seq-id {identity} --cov-mode {cov_mode} --cluster-mode 2 --cluster-reassign -v 3 --threads {self.threads}"
        )
        self._run_command(command, log_file)

    def createtsv(self, tsv_file='result_seq_clu.tsv', log_file='createTsv.log'):
        command = f"{self.mmseqs2_path} createtsv {self.dir}/db {self.dir}/db {self.dir}/db_clu {self.dir}/{tsv_file}"
        self._run_command(command, log_file)
        if self.cleanup:
            print("Cleaning up intermediate files")
            self._cleanup_intermediate_files(tsv_file)

    def _cleanup_intermediate_files(self, tsv_file):
        files_to_keep = {self.dir / tsv_file}
        for file in self.dir.iterdir():
            if file not in files_to_keep:
                try:
                    if file.is_file():
                        file.unlink()
                    elif file.is_dir():
                        shutil.rmtree(file)
                except Exception as e:
                    print(f"Error deleting file {file}: {e}")

    def fasta2representativeseq(self, fasta_file : Path, cov : float, iden : float, cov_mode : int = 0) -> dict:
        
        self.createdb(fasta_file)
        self.cluster(coverage=cov, identity=iden, cov_mode=cov_mode)
        self.createtsv()

        with open(self.dir / "result_seq_clu.tsv") as tsv:
            representatives = {line.split()[0] for line in tsv.readlines()}
            
        full_fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        representatives_fasta = {k: full_fasta[k] for k in representatives}
        
        return representatives_fasta
    
    def fasta2representative_mobidb(self, fasta_file : Path, cov : float, iden : float, cov_mode : int = 0) -> dict:

        self.createdb(fasta_file)
        self.cluster(coverage=cov, identity=iden, cov_mode=cov_mode)
        self.createtsv()

        priority = {"curated": 3, "derived": 2, "homology": 1, "predicted" : 0}

        clusters = {}
        with open(self.dir / "result_seq_clu.tsv") as tsv:
            
            for line in tsv.readlines():
                representative, member = line.split()
                if representative not in clusters:
                    clusters[representative] = list()
                clusters[representative].append(member)
                clusters[representative].append(representative) # representative is also a member of its own cluster 

        representatives = list() # Not a set, we expect MMseqs to only allocate each sequence once in one cluster
        for _, members in clusters.items():
            if len(members) == 1:
                representatives.append(members[0])
            else:
                best_member = members[0]
                for member in members[1:]:
                    if priority[member.split("|")[1]] > priority[best_member.split("|")[1]]:
                        best_member = member
                    elif priority[member.split("|")[1]] == priority[best_member.split("|")[1]]:
                        member_len = int(member.split("|")[2].split("-")[1]) - int(member.split("|")[2].split("-")[0])
                        best_member_len = int(best_member.split("|")[2].split("-")[1]) - int(best_member.split("|")[2].split("-")[0])
                        if member_len > best_member_len:
                            best_member = member
                representatives.append(best_member)


        full_fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        representatives_fasta = {k: full_fasta[k] for k in representatives}

        return representatives_fasta


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--fasta', help='Path to the input FASTA file')
    parser.add_argument('--writing_dir', default = ".", help='Directory to write the representative sequences')
    parser.add_argument('--cov', type=float, default = 0.7, help='Coverage threshold for clustering')
    parser.add_argument('--iden', type=float, default = 0.3, help='Identity threshold for clustering')
    parser.add_argument('--cov_mode', type=int, default=0, help='Coverage mode for clustering')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use')
    parser.add_argument('--cleanup', action='store_true', default=True, help='Cleanup intermediate files')
    parser.add_argument('--mmseqs_path', default="mmseqs", help='Path to the MMseqs2 binary')

    args = parser.parse_args()

    mmseqs2_api = MMseqs2(threads=args.threads, cleanup=args.cleanup, mmseqs2_path=args.mmseqs_path)
    
    representatives = mmseqs2_api.fasta2representativeseq(args.fasta, args.cov, args.iden, args.cov_mode)

    SeqIO.write(representatives.values(), f"{args.writing_dir}/representative.fasta", "fasta")
    
    

   
