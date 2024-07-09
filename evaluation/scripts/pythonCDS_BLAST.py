from Bio import SeqIO
import gzip
import sys

def makecds(inpath,outprefix):


    with gzip.open(inpath,"rt") as gbff:
        with open(outprefix+"_genome.fasta","w") as genome_out:
            with open( outprefix+"_cds.bed","w") as bed_out:

                for entry in SeqIO.parse(gbff, "genbank"):
                    sequenceid = entry.id
                    SeqIO.write(entry,genome_out,"fasta")
                    cds_blocks=[]
                    for feature in entry.features:
                        if feature.type == "CDS":
                            for seqElem in feature.location.parts:
                                cds_blocks.append((sequenceid,seqElem.nofuzzy_start,seqElem.nofuzzy_end))
                    for sequenceid,seqstart,seqend in sortmergeBed(cds_blocks):
                        bed_out.write(f"{sequenceid}\t{seqstart}\t{seqend}\n")

def sortmergeBed(bedtuplelist):
    bedtuplelist.sort()
    if len(bedtuplelist)<2:
        return(bedtuplelist)
    r = []
    for tuple in bedtuplelist:
        assert tuple[2]>=tuple[1]
        if not r:
            r.append(tuple)
        else:
            lastelem = r[-1]
            if lastelem[2]>=tuple[1] and tuple[0]==lastelem[0]:
                new = (lastelem[0],lastelem[1],max(lastelem[2],tuple[2]))
              #  print("merged", "%s,%d,%d" % lastelem, "with", "%s%d%d" % tuple, "to", "%s%d%d" % new)
                r[-1]=new
               # print("---", r[-1])
            else:
                r.append(tuple)
    return(r)


def blastoutfmt6_to_bed3(genomebed,tblastnresult,proteinsequences, outfile):
    genomecds = []
    with open(genomebed,"r") as bed:
        for elem in bed:
            if elem:
                sid, start,stop = elem.strip().split("\t")
                genomecds.append((sid,int(start),int(stop)))
    # Protein Sequences
    protseqlengths = load_ProtSequences(proteinsequences)
    
    
    # Define Columns in Blast OUTFMT6 format
    refstartpos = 8
    refstoppos = 9
    seqidentitypos = 2
    querystartpos = 6
    querystoppos = 7
    
    # Define Thresholds
    percentIdentity = 90
    lengthIdentity = 0.75

    with open(tblastnresult,"r") as blasttsv:
        for hit in blasttsv:
            hit_split = hit.strip().split()
            if len(hit_split):
                sid = hit_split[1]
                queryid = hit_split[0]
                
                if not queryid in protseqlengths:
                    print (queryid, "not found in", proteinsequences )
                required_protlength = protseqlengths[queryid]*lengthIdentity
                hitlength = abs(int(hit_split[querystoppos])-int(hit_split[querystartpos]))
                
                pident = float (hit_split[seqidentitypos])
                if pident >= percentIdentity and hitlength >= required_protlength:
                    start = int(hit_split[8])
                    stop = int(hit_split[9])


                    minstart = min(start,stop)-1 # to getBED 0 based
                    maxstop = max(start,stop)   # no substraction since we need end+1 

                    genomecds.append((sid,minstart,maxstop))
    
    with open(outfile,"w") as out:
        for elem in sortmergeBed(genomecds):
            out.write(f"{elem[0]}\t{elem[1]}\t{elem[2]}\n" )




#############################################
# Load Protein Sequences and Store lengths
#############################################

def load_ProtSequences(fastafile,name=""):
    with open (fastafile,"r") as proteinsequences:
        seqname="null"
        proteinlengths = dict()
        for line in proteinsequences:
            # if new entry begin
            if line.startswith(">"):
                seqname = line.lstrip(">").split(" ")[0].strip()
                
                if seqname in proteinlengths:
                    print(f"Warning: {seqname} multiple times in {fastafile}")
                else:
                    proteinlengths[seqname]=0
            else:
                proteinlengths[seqname] += len(line.strip())
    print(f"{name+' ' if name else ''} Loaded {len(proteinlengths)} unique Protein lengths")
    return(proteinlengths)










if sys.argv[1]== "makecds":
    print(sys.argv)
    makecds(sys.argv[2],sys.argv[3])
if sys.argv[1]=="blast2bed":
    print("blast2bed")
    blastoutfmt6_to_bed3(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
