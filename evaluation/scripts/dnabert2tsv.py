import pandas
import numpy as np
import datetime

def dnabert2tsv(classinput,numpyresult,output, logfile = sys.stdout):
    print("Date:\t\t\t",datetime.datetime.now().isoformat(), file =logfile)
    print("Groundtruth:\t\t",classinput, file =logfile)
    print("DNAbert Prediction:\t",numpyresult, file=logfile)
    print("Anotated Result:\t", output, file = logfile)
    with open(classinput,"r") as template:
        with open(output,"w") as opf:
            header = template.readline().strip().split()[1:]
            header.append("Prediction")
            headerline = "\t".join(header) + "\n"
            opf.write(headerline)
            pred = np.load(numpyresult)
        
            for index, line in enumerate(template):
                tl = line.strip().split("\t")[1:]

                tl.append(str(pred[index]))
                tls = "\t".join(tl) + "\n"
                opf.write(tls)

with open(snakemake.log[0],"w") as log:
    print(snakemake.wildcards["sample"],snakemake.input.keys(),file=log)

    dnabert2tsv(snakemake.input["groundtruth"],snakemake.input["predictionfile"],snakemake.output[0],log)
