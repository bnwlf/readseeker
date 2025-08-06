import os, csv, glob

header = ["Sequence","Label","Name","MappingQuality","Referenz","Sample","ReadNotInBlastHits"]
coding = 1
noncoding = 0


outputfile = "../benchmark/human_bm/demo_coding_vs_intergenomic_seqs/dev.tsv"
ofh = open(outputfile,"w")
ofh_csv = csv.writer(ofh,delimiter='\t')
ofh_csv.writerow(header)

fastafh = open("../benchmark/human_bm/demo_coding_vs_intergenomic_seqs/demo_coding_vs_intergenomic_seqs_reads.fasta","w")


def file2entry(file,label,fa):
    with open(file,"r") as seq:
        sequence = seq.readlines(1)[0].strip()
        kmer_seq = [sequence[i:i+6] for i in range(len(sequence)-5)]
        f = "/".join(file.split("/")[-4:])
        fastafh.write(f">{f}-{label}\n{sequence}\n")
    return([" ".join(kmer_seq),label,file,0,"human_bm","demo_coding_vs_intergenomic_seqs_test",0])


for file in sorted(glob.glob(os.path.expanduser("~/.genomic_benchmarks/demo_coding_vs_intergenomic_seqs/test/coding_seqs/*.txt")), key = lambda fn: int(os.path.basename(fn)[:-4])):
	ofh_csv.writerow(file2entry(file,coding,fastafh))

for file in sorted(glob.glob(os.path.expanduser("~/.genomic_benchmarks/demo_coding_vs_intergenomic_seqs/test/intergenomic_seqs/*.txt")), key = lambda fn: int(os.path.basename(fn)[:-4])):
    ofh_csv.writerow(file2entry(file,noncoding,fastafh))

del ofh_csv
ofh.close()
fastafh.close()
