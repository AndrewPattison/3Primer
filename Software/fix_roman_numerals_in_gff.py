#!/usr/bin/python
import sys
import roman
# you might have to pip install roman if it doesn't work right away

# Make any roman numerals into ints

def fix_gff_file(gff_file):
    with open("fixed_gff_file.gff", 'w') as out_f:
        with open(gff_file, 'r') as in_f:
            for line in in_f:
                if line[0] == "#":
                    continue
                else:
                    splut = line.split("\t")
                    try:
                        splut[0] = roman.fromRoman(splut[0])
                    except:
                        pass
                    out_f.write(str(splut[0]) + "\t" + str(splut[1]) + "\t" + str(splut[2]) + "\t" + splut[3] + 
                        "\t" + splut[4] + "\t" + splut[5] + "\t" + splut[6] + "\t" + splut[7] + "\t" + splut[8]) 

fix_gff_file(sys.argv[1])