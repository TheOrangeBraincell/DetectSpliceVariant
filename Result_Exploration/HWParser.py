# HW Parser
# 19.02.24
# Mirjam Müller

import sys
import re

input_file = sys.argv[1]
current_variant = ""
pvalue_found = False
results_dict = dict()
with open(input_file, "r") as infile, open("HW_results_parsed.tsv", "w") as out:
    for line in infile:
        # skip header
        if line.startswith("l"):
            out.write("Location\tp-value\n")
            continue
        variant = line.split("\t")[0]
        if len(line.split("\t")) > 1:
            infostring = line.strip("\n").split("\t")[1]
        else:
            # no info string on this line.
            continue
        if variant != current_variant:
            # we found a new entry
            # if we did not find a p-value for the previous variant, write it to console.
            if pvalue_found == False and current_variant != "":
                print("No p value found!")
                print(current_variant)
            if current_variant != "":
                #write output
                out.write(current_variant+"\t"+results_dict[current_variant]+"\n")
            # reset variables
            current_variant=variant
            pvalue_found = False

        if "p-value" in infostring:
            # found p-value, use regex to extract it.
            pvalue = re.search(r"p-value [=<>] (\d+.?\d*e?-?\d*)", infostring)[0].split(" ")[2]
            pvalue_found = True
            if current_variant not in results_dict:
                results_dict[current_variant]=pvalue
            else:
                #this should not happen. we dont want to find several pvalues
                print("Several p-values found?")
                print(current_variant)
                print(results_dict[current_variant])
                print(infostring)
        elif "No variant alleles observed" in infostring:
            # There has been no test possible. So the p-value is a NAN.
            pvalue_found = True
            if current_variant not in results_dict:
                results_dict[current_variant] = "NAN"
            else:
                # this should not happen. we dont want to find several pvalues
                print("p-values and a failed test?")
                print(current_variant)
                print(results_dict[current_variant])
                print(infostring)

    #Write last line
    if pvalue_found == False and current_variant != "":
        print("No p value found!")
        print(current_variant)
    if current_variant != "":
        # write output
        out.write(current_variant + "\t" + results_dict[current_variant]+"\n")