import sys
import os

input_files = os.listdir(sys.argv[1])

outfile = open(sys.argv[2], "w")
outfile.write("filetype\tEvent\tcount\n")
for file in input_files:
    with open(sys.argv[1]+file, "r") as f:
        for line in f:
            if line.startswith("Event"):
                filetype="PSI"
                continue
            elif line.startswith("Location"):
                filetype="Genotype"
                continue

            entries = line.strip().split("\t")
            count = 0
            for e in entries[1:]:
                if filetype=="PSI":
                    if e!="NAN":
                        count+=1
                elif filetype=="Genotype":
                    if e!="ND" and e!="NE":
                        count+=1
                else:
                    print("Error: Unknown filetype")
                    sys.exit(1)
            
            outfile.write(filetype + "\t" + entries[0] + "\t" + str(count) + "\n")

outfile.close()

