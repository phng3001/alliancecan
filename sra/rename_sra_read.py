# P=NP

import sys

def rename_reads(sra_prefix, input_fastq, output_fastq):
    with open(input_fastq, 'r') as fastq, open(output_fastq, 'w') as output:
        for line in fastq:
            if line.startswith(f"@{sra_prefix}"):
                tokens = line.split(" ")
                readID = tokens[0][0:-2]  # Remove the last .1 or .2
                newline = readID + " " + " ".join(tokens[1:])
                output.write(newline)
            else:
                output.write(line)

import sys

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python rename_read.py <sra_prefix> <input_fastq> <output_fastq>")
        sys.exit(1)
    
    sra_prefix = sys.argv[1]
    input_fastq = sys.argv[2]
    output_fastq = sys.argv[3]
    
    rename_reads(sra_prefix, input_fastq, output_fastq)



#input
#@ERR433977.1.1 1 length=100
#ACCATAGATAGCAAGAATTTTCCACAAGCTGTGAAAAATAATCCACAATGTTACCAAACTTTATCCACAGGTTGGGGATAAAAGAAGAAATTATTGATTT
#+ERR433977.1.1 1 length=100
#?<??=DD>4CF>DG9ADE9E:?EHG<F8@?AA?<C@??)9:?F9*?FF@EB999BFGG8BB8B==8=7=8=.@35'5;?@:@>3;@;5?5;(:>;@5535

#@ERR433977.1.2 1 length=100
#AAGGTTATCCACTATGTTTTTCGATAAAAAGCTTAATAAATCAATAATTTCTTCTTTTATCCCCAACCTGTGGATAAAGTTTGGTAACATTGTGGATTAT
#+ERR433977.1.2 1 length=100
#=<?D4=AD:DDD>E@ECEECDCDAFEFIEEI;<EDE<CDDF<BDDADEIDABDDCD>DDEC4BDDICDEIEDAA?CDDD.);7;66.;6>>5>A?AA>>>



#output
#@ERR433977.1 1 length=100
#ACCATAGATAGCAAGAATTTTCCACAAGCTGTGAAAAATAATCCACAATGTTACCAAACTTTATCCACAGGTTGGGGATAAAAGAAGAAATTATTGATTT
#+ERR433977.1.1 1 length=100
#?<??=DD>4CF>DG9ADE9E:?EHG<F8@?AA?<C@??)9:?F9*?FF@EB999BFGG8BB8B==8=7=8=.@35'5;?@:@>3;@;5?5;(:>;@5535

#@ERR433977.1 1 length=100
#AAGGTTATCCACTATGTTTTTCGATAAAAAGCTTAATAAATCAATAATTTCTTCTTTTATCCCCAACCTGTGGATAAAGTTTGGTAACATTGTGGATTAT
#+ERR433977.1.2 1 length=100
#=<?D4=AD:DDD>E@ECEECDCDAFEFIEEI;<EDE<CDDF<BDDADEIDABDDCD>DDEC4BDDICDEIEDAA?CDDD.);7;66.;6>>5>A?AA>>>
