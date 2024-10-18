# Z-NOME_HUNT

## What is this?
This github respository contains source code and a binary for a modernized version of the algorithim Z-Hunt, now able to perform analysis of large genomic files in a user-friendly manner. It returns a "Z-Score" (not related to the common statistical Z-Score) for each nucleotide in a given sequence to rate its propensity to convert to Z-Form, which is then paired with the location of the nucleotide in the genome and outputted as a BedGraph. 
Z-NOME_HUNT runs as a command line interface tool. It takes four necessary arguments:

`window_size` : The size of the analysis window after a nucleotide. Larger values will take longer to run.
`min`: Minimum window size acceptable for a Z-forming region.

`max`: The maximum acceptable size for a z forming region. Larger deltas between minimum and maximum will increase the runtime.

`file-name`: Your file of interest. See formatting notes.

and two optional arguments:

`chromosome_name`: The chromosome number or name your region of interest comes from.

`start_site`: The starting point for your sequence.

The original Z-HUNT and MZ-HUNT code that formed the basis for this project can be found here:

https://github.com/Ho-Lab-Colostate/zhunt/tree/master

Useful background reading can be found here:

https://github.com/Ho-Lab-Colostate/zhunt/tree/master/docs

https://link.springer.com/protocol/10.1007/978-1-0716-3084-6_14

## Quick start
Download z-nome_hunt.c in zhunt_genome/source, then (if using gcc):
`gcc -o z-nome_hunt z-nome_hunt.c -lm`

then
`./z-nome_hunt 12 6 12 <your_file>`

A new file will be created, <yourfile>.Z-SCORE.bedgraph. This can then be converted into a BigWig file by using `bedgraphToBigWig` or `wigToBigWig`.

By default, Z-NOME_HUNT assigns your filename (with .txt, .fa, or .fna extension removed) as the chromosome name in the output bedgraph file, and initializes start site as 0. If you want to run the algorithm with a different start site, use

`./z-nome_hunt 12 6 12 <your_file> <chromosome_name> <start_site>`


## Important formatting notes for input file
Your input file should 
# Contain no comments. 
The only information in the file should be its sequence. Downloaded chromosomes will usually have a first line with > then information. You can run

`sed -i '1d' <yourfile>`

in bash to remove this first line.

# Be a .txt, .fa, or .fna file, or have its file extension removed. 
Formats besides these extensions or compressed sequence files may not work properly.


'N's in the genome and bases neighboring 'N's are assigned an arbitrary very low Z-Score. Z-scores for bases within 6-12 bases of an 'N' are therefore not trustworthy and should be ignored.
