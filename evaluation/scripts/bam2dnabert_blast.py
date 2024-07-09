import argparse, sys, random, pysam
import numpy as np
import datetime


def bam2props(filename,minimal_readlength=250):
    readpropsdict = dict()
    for read in  pysam.AlignmentFile(filename, "rb"):
        readname = "{}_{}".format(read.qname, read.is_read2+1)
        if len(read.seq) >= minimal_readlength:
            if readname not in readpropsdict:
                readpropsdict[readname] = -1
            readpropsdict[readname] = max(readpropsdict[readname],read.mapping_quality)
    return readpropsdict

def filtered2tsvrows(pullnames, filename, typ, k=6,reference=None,sample=None,blastnames=set()):
    rows = []
    assert typ in ["cds","noncds"]
    label ={"cds":1,"noncds":0}.get(typ)
    
    st = "" if not sample else f"\t{sample}"
    rt = "" if not reference else f"\t{reference}"
    
    for read in  pysam.AlignmentFile(filename, "rb"):
        readname = "{}_{}".format(read.qname, read.is_read2+1)
        if readname in pullnames:
            sequence = ' '.join([read.query[x:x+k] for x in range(len(read.query)-k+1)])
            rows.append(f"{sequence}\t{label}\t{readname}\t{read.mapping_quality}{rt}{st}\t{1 if readname in blastnames else 0}\n")
            pullnames.remove(readname) # Avoid selecting read twice
    return(rows)  

def bams2dnabert(cds,noncds,outputfile,qualitycutoff = 40, balance = False, balanceseed = 42,logfile = sys.stdout, referencename=None, samplename=None,minimalreadlength=140, blastReads = None):
    print("CDS File:\t\t",cds, file = logfile)
    print("nonCDS mapping File:\t\t", noncds,file = logfile)
    print("nonCDS mappingblast File",blastReads if blastReads else "-",file = logfile)
    print("Outputfile:\t\t", outputfile, file = logfile)
    print("Quality Cutoff:\t\t", qualitycutoff, file = logfile)
    print("CDS/NONCDS balanced:\t", balance, file = logfile)
    print("balanceseed:\t\t", balanceseed, file = logfile)
    if (referencename):
        print("Reference:\t\t", referencename, file = logfile)
    if (samplename):
        print("Sample:\t\t\t", samplename, file = logfile)
    print("Date:\t\t\t", datetime.datetime.now().isoformat(), file = logfile)
    
    #Store Readnames and Mappingqualities per Read
    cdsreadsProps = bam2props(cds,minimalreadlength)
    noncdsreadsProps = bam2props(noncds,minimalreadlength)
    if blastReads:
        blastReadsProps = bam2props(blastReads)

    # Exclude Reads occuring in both categories
    uniq_cds_Props = {readname: quality for readname,quality in cdsreadsProps.items() if readname not in noncdsreadsProps}
    uniq_noncds_Props = {readname: quality for readname,quality in noncdsreadsProps.items() if readname not in cdsreadsProps}
    intersection = set(cdsreadsProps.keys()).intersection(noncdsreadsProps)
    print("uniquecdsProps", len(uniq_cds_Props), file = logfile)
    print("unique NON cdsProps", len(uniq_noncds_Props), file = logfile)
    print("intersection",len(intersection), file = logfile)

    if blastReads:
        print("uniq Blastreads",len(blastReadsProps),file = logfile)
    
    
    hqreads_cds = {readname:quality for readname, quality in uniq_cds_Props.items() if quality >= qualitycutoff }
    hqreads_noncds = {readname:quality for readname, quality in uniq_noncds_Props.items() if quality >= qualitycutoff }
    if blastReads:
        hqreads_blastnoncds =  {readname:quality for readname, quality in  blastReadsProps.items() if quality >= qualitycutoff }
    else:
        hqreads_blastnoncds = dict()
    print("\n", file = logfile)
    print("_________Quality Filtered Reads_________", file = logfile)
    print("QF CDS:", len(hqreads_cds), file = logfile)
    print("QF nonCDS:", len(hqreads_noncds), file = logfile)
    if blastReads:
        print("QF nonCDS blast:", len(hqreads_blastnoncds), file = logfile)
    
    if balance:
        maxreads = min(len(hqreads_cds),len(hqreads_noncds))
        

        if (len(hqreads_cds) and len(hqreads_noncds)):
            print("""
_________Balancing Losses__________
Balanced CDS {0}, {1}% loss
Balanced nonCDS {0}, {2}% loss""".format(
            maxreads,
            100-((maxreads*100)/len(hqreads_cds)),
            100-((maxreads*100)/len(hqreads_noncds))
        ), file = logfile)
        random.seed(balanceseed)
        
        # Sample random readnames
        hqreadnames_sorted_cds = set(random.sample(sorted(hqreads_cds.keys()),maxreads))
        hqreadnames_sorted_noncds = set(random.sample(sorted(hqreads_noncds.keys()),maxreads))

        
    
        with open(outputfile, "w") as out:
            rt = "" if not referencename else f"\tReferenz"
            st = "" if not samplename else f"\tSample"
            out.write(f"Sequence\tLabel\tName\tMappingQuality{rt}{st}\tReadNotInBlastHits\n")
            out.writelines(filtered2tsvrows(
                hqreadnames_sorted_cds,
                cds,
                "cds",
                6,
                referencename,
                samplename),)
            
            out.writelines(filtered2tsvrows(
                hqreadnames_sorted_noncds,
                noncds,
                "noncds",
                6,
                referencename,
                samplename))
    else:
      with open(outputfile, "w") as out:
            rt = "" if not referencename else f"\tReferenz"
            st = "" if not samplename else f"\tSample"
            out.write(f"Sequence\tLabel\tName\tMappingQuality{rt}{st}\n")
            out.writelines(filtered2tsvrows(
                set(hqreads_cds),
                cds,
                "cds",
                6,
                referencename,
                samplename))
            
            out.writelines(filtered2tsvrows(
                set(hqreads_noncds),
                noncds,
                "noncds",
                6,
                referencename,
                samplename,set(hqreads_blastnoncds.keys())))
    print(f"_________Done_________\n{datetime.datetime.now().isoformat()}", file = logfile)
        
    
with open(snakemake.log[0], "w") as log:
  log.write("Wildcards:\n")
  for key,value in snakemake.wildcards.items():
      log.write(f"\t{key}:\t\t{value}\n")
      log.write("______________________\n")
  blastreads = None
  if "mapblastnonCDS" in snakemake.input:
      blastreads = snakemake.input["mapblastnonCDS"]
  bams2dnabert(snakemake.input["cds"],
             snakemake.input["mapnonCDS"],
             snakemake.output[0],
             snakemake.params["minimalmappingquality"],
             snakemake.params["balanced"],
             referencename=snakemake.wildcards["reference"],
             samplename=snakemake.wildcards["sample"],
             blastReads=blastreads,
             logfile=log
		)
