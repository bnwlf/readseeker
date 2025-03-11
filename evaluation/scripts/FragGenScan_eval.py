import sys


def load_original_fasta(fastafile):
    with open(fastafile) as f:
        readdict = dict()
        name = None
        lenlength = 0
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    readdict[name] = lenlength
                name = line[1:].strip()
                lenlength = 0
            else:
                lenlength += len(line.strip())
        if name is not None:
            readdict[name] = lenlength
        return readdict

def summcds_regions(readlendict, fastaresult,outfile = sys.stdout, samplename = "unknown",reference="unknown"):

    cdssums = {key: 0 for key in readlendict.keys()}

    with open(fastaresult, "r") as f:
        for line in f:


            if line.startswith(">"):
                header_elements = line[1:].strip().split("_")

                first = int(header_elements[-3])-1
                last = int(header_elements[-2])
                name = "_".join(header_elements[:-3])
                cdssums[name] = abs(first - last)


    with open(outfile,"w") as f:
        f.write("Reference\tSample\tName\tcdsbases\ttotallength\tfraction\tPredictedLabel\tLabel\n")
        template ="%s\t%s\t%s\t%d\t%d\t%f\t%d\t%s\n"
        for readname, length in readlendict.items():

            fraction = float(cdssums.get(readname,0))/length
            classification = 0 if fraction < 0.5 else 1

            ref_class = readname.split("-")[-1]

            f.write(template % (reference, samplename, readname, cdssums.get(readname), length,fraction,classification, ref_class))


rld = load_original_fasta(snakemake.input["fasta"])
summcds_regions(rld, snakemake.input["ffn"], snakemake.output[0], samplename=snakemake.wildcards["sample"],reference=snakemake.wildcards["reference"])






