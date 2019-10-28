# bioinformatics-one-liners
my collection of bioinformatics one liners that is useful in my day-to-day work

### I came across the bioinformatics one-liners on the [biostar](https://www.biostars.org/p/142545/) forum and gathered them here.
I also added some of my own tricks

05/21/2015.



####  get the sequences length distribution form a fastq file using awk

```bash
zcat file.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'  
```
#### Reverse complement a sequence (I use that a lot when I need to design primers)

```
echo 'ATTGCTATGCTNNNT' | rev | tr 'ACTG' 'TGAC'
```

#### split a multifasta file into single ones with csplit:

```bash
csplit -z -q -n 4 -f sequence_ sequences.fasta /\>/ {*}  
```
#### Split a multi-FASTA file into individual FASTA files by awk

```bash
awk '/^>/{s=++d".fa"} {print > s}' multi.fa
```

#### linearize multiline fasta

```bash
cat file.fasta | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}'
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' file.fa
```
#### fastq2fasta

```bash
zcat file.fastq.gz | paste - - - - | perl -ane 'print ">$F[0]\n$F[2]\n";' | gzip -c > file.fasta.gz
```
####  bam2bed

```bash
samtools view file.bam | perl -F'\t' -ane '$strand=($F[1]&16)?"-":"+";$length=1;$tmp=$F[5];$tmp =~ s/(\d+)[MD]/$length+=$1/eg;print "$F[2]\t$F[3]\t".($F[3]+$length)."\t$F[0]\t0\t$strand\n";' > file.bed
```

#### bam2wig

```bash
samtools mpileup -BQ0 file.sorted.bam | perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print "fixedStep chrom=$c start=$start step=1 span=1\n";}$_ = $depth."\n";($lastC, $lastStart) = ($c, $start);' | gzip -c > file.wig.gz
```

#### Number of reads in a fastq file

```bash
cat file.fq | echo $((`wc -l`/4))
```
#### Single line fasta file to multi-line fasta of 60 characteres each line

```bash
awk -v FS= '/^>/{print;next}{for (i=0;i<=NF/60;i++) {for (j=1;j<=60;j++) printf "%s", $(i*60 +j); print ""}}' file

fold -w 60 file
```

#### Sequence length of every entry in a multifasta file

```bash
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' file.fa
```
#### Reproducible subsampling of a FASTQ file. srand() is the seed for the random number generator - keeps the subsampling the same when the script is run multiple times.  0.01 is the % of reads to output.

```bash
cat file.fq | paste - - - - | awk 'BEGIN{srand(1234)}{if(rand() < 0.01) print $0}' | tr '\t' '\n' > out.fq
```
#### or look at the Hengli's Seqtk 

#### Deinterleaving a FASTQ:

```bash
cat file.fq | paste - - - - - - - - | tee >(cut -f1-4 | tr '\t'  
'\n' > out1.fq) | cut -f5-8 | tr '\t' '\n' > out2.fq
```

#### Using mpileup for a whole genome can take forever. So, handling each chromosome separately and parallely running them on several cores will speed up your pipeline. Using xargs you can easily realize it.  
#### Example usage of xargs (-P is the number of parallel processes started - don't use more than the number of cores you have available):

```basg
samtools view -H yourFile.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d 100000 -uf yourGenome.fa -r {} yourFile.bam | bcftools view -vcg - > tmp.{}.vcf"
```

#### To merge the results afterwards, you might want to do something like this:

```bash
samtools view -H yourFile.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | perl -ane 'system("cat tmp.$F[0].bcf >> yourFile.vcf");'
```

#### split large file by id/label/column

```bash
awk '{print >> $1; close($1)}' input_file
```
#### split a bed file by chromosome:

```bash
cat nexterarapidcapture_exome_targetedregions_v1.2.bed | sort -k1,1 -k2,2n | sed 's/^chr//' | awk '{close(f);f=$1}{print > f".bed"}'

#or
awk '{print $0 >> $1".bed"}' example.bed
```

#### sort vcf file with header

```bash
cat my.vcf | awk '$0~"^#" { print $0; next } { print $0 | "sort -k1,1V -k2,2n" }'
```
#### Rename a file, bash string manipulation

```bash
for file in *gz
do zcat $file > ${file/bed.gz/bed}
```

#### gnu sed print invisible characters

```bash
cat my_file | sed -n 'l'
cat -A
```

#### exit a dead ssh session
`~.`

#### copy large files, copy the from_dir directory inside the to_dir directory

```bash
rsync -av from_dir  to_dir

## copy every file inside the frm_dir to to_dir
rsync -av from_dir/ to_dir

##re-copy the files avoiding completed ones:

rsync -avhP /from/dir /to/dir
```

#### make directory using the current date

```bash
mkdir $(date +%F)
```
#### all the folders' size in the current folder (GNU du)

```bash
du -h --max-depth=1
```

### this one is a bit different, try it and see the difference
`du -ch`

#### the total size of current directory
`du -sh .`

#### disk usage
`df -h`

#### the column names of the file, install csvkit https://csvkit.readthedocs.org/en/0.9.1/
`csvcut -n`

#### open top with human readable size in Mb, Gb. install htop for better visualization
`top -M`

#### how many memeory are used in Gb
`free -mg`

#### print out unique rows based on the first and second column
`awk '!a[$1,$2]++' input_file`

`sort -u -k1,2 file`
It will sort based on unique first and second column

#### do not wrap the lines using less
`less -S`

#### pretty output
```bash
fold -w 60
cat file.txt | column -t | less -S
```
#### pass tab as delimiter http://unix.stackexchange.com/questions/46910/is-it-a-bug-for-join-with-t-t
`-t $'\t'`

#### awk with the first line printed always
`awk ' NR ==1 || ($10 > 1 && $11 > 0 && $18 > 0.001)'  input_file`

#### delete blank lines with sed
`sed /^$/d`

#### delete the last line
`sed $d`

awk to join files based on several columns

my [github repo](https://github.com/crazyhottommy/scripts-general-use/blob/master/Shell/Awk_anotates_vcf_with_bed.ipynb)

```
### select lines from a file based on columns in another file
## http://unix.stackexchange.com/questions/134829/compare-two-columns-of-different-files-and-print-if-it-matches
awk -F"\t" 'NR==FNR{a[$1$2$3]++;next};a[$1$2$3] > 0' file2 file1 

```

Finally learned about the !$ in unix: take the last thing (word) from the previous command.   
`echo hello, world; echo !$` gives 'world'


Create a script of the last executed command:  
`echo "!!" > foo.sh`

Reuse all parameter of the previous command line:  
`!*`

find bam in current folder (search recursively) and copy it to a new directory using 5 CPUs    
`find . -name "*bam" | xargs -P5 -I{} rsync -av {} dest_dir`

`ls -X`  will group files by extension.

loop through all the chromosomes

```bash
for i in {1..22} X Y 
do
  echo $i
done
```

for i in in `{01..22}` will expand to 01 02 ...


change every other newline to tab:

`paste` is used to concatenate corresponding lines from files: paste file1 file2 file3 .... If one of the "file" arguments is "-", then lines are read from standard input. If there are 2 "-" arguments, then paste takes 2 lines from stdin. And so on.

```bash
cat test.txt  
0    ATTTTATTNGAAATAGTAGTGGG
0    CTCCCAAAATACTAAAATTATAA
1    TTTTAGTTATTTANGAGGTTGAG
1    CNTAATCTTAACTCACTACAACC
2    TTATAATTTTAGTATTTTGGGAG
2    CATATTAACCAAACTAATCTTAA
3    GGTTAATATGGTGAAATTTAAT
3    ACCTCAACCTCNTAAATAACTAA

cat test.txt| paste - -                               
0    ATTTTATTNGAAATAGTAGTGGG    0    CTCCCAAAATACTAAAATTATAA
1    TTTTAGTTATTTANGAGGTTGAG    1    CNTAATCTTAACTCACTACAACC
2    TTATAATTTTAGTATTTTGGGAG    2    CATATTAACCAAACTAATCTTAA
3    GGTTAATATGGTGAAATTTAAT     3    ACCTCAACCTCNTAAATAACTAA
```

ORS: output record seperator in `awk`
`var=condition?condition_if_true:condition_if_false is the ternary operator.`

```bash
cat test.txt| awk 'ORS=NR%2?"\t":"\n"'          

0    ATTTTATTNGAAATAGTAGTGGG    0    CTCCCAAAATACTAAAATTATAA
1    TTTTAGTTATTTANGAGGTTGAG    1    CNTAATCTTAACTCACTACAACC
2    TTATAATTTTAGTATTTTGGGAG    2    CATATTAACCAAACTAATCTTAA
3    GGTTAATATGGTGAAATTTAAT     3    ACCTCAACCTCNTAAATAACTAA

```

#### awk
We can also use the concept of a conditional operator in print statement of the form print CONDITION ? PRINT_IF_TRUE_TEXT : PRINT_IF_FALSE_TEXT. For example, in the code below, we identify sequences with lengths > 14:

```bash
cat data/test.tsv
blah_C1	ACTGTCTGTCACTGTGTTGTGATGTTGTGTGTG
blah_C2	ACTTTATATATT
blah_C3	ACTTATATATATATA
blah_C4	ACTTATATATATATA
blah_C5	ACTTTATATATT	

awk '{print (length($2)>14) ? $0">14" : $0"<=14";}' data/test.tsv
blah_C1	ACTGTCTGTCACTGTGTTGTGATGTTGTGTGTG>14
blah_C2	ACTTTATATATT<=14
blah_C3	ACTTATATATATATA>14
blah_C4	ACTTATATATATATA>14
blah_C5	ACTTTATATATT<=14

awk 'NR==3{print "";next}{printf $1"\t"}{print $1}' data/test.tsv
blah_C1	blah_C1
blah_C2	blah_C2

blah_C4	blah_C4
blah_C5	blah_C5

```
You can also use getline to load the contents of another file in addition to the one you are reading, for example, in the statement given below, the while loop will load each line from test.tsv into k until no more lines are to be read:
```bash
awk 'BEGIN{while((getline k <"data/test.tsv")>0) print "BEGIN:"k}{print}' data/test.tsv
BEGIN:blah_C1	ACTGTCTGTCACTGTGTTGTGATGTTGTGTGTG
BEGIN:blah_C2	ACTTTATATATT
BEGIN:blah_C3	ACTTATATATATATA
BEGIN:blah_C4	ACTTATATATATATA
BEGIN:blah_C5	ACTTTATATATT
blah_C1	ACTGTCTGTCACTGTGTTGTGATGTTGTGTGTG
blah_C2	ACTTTATATATT
blah_C3	ACTTATATATATATA
blah_C4	ACTTATATATATATA
blah_C5	ACTTTATATATT
```
#### merge multiple fasta sequences in two files into a single file line by line
see [post](https://www.biostars.org/p/204336/#204380)  

`linearize.awk:`  

```bash
/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}
```

```bash
paste <(awk -f linearize.awk file1.fa ) <(awk -f linearize.awk file2.fa  )| tr "\t" "\n"
```

#### grep fastq reads containing a pattern but maintain the fastq format

```bash
grep -A 2 -B 1 'AAGTTGATAACGGACTAGCCTTATTTT' file.fq | sed '/^--$/d' > out.fq

# or
zcat reads.fq.gz \
| paste - - - - \
| awk -v FS="\t" -v OFS="\n" '$2 ~ "AAGTTGATAACGGACTAGCCTTATTTT" {print $1, $2, $3, $4}' \
| gzip > filtered.fq.gz
```

#### count how many columns of a tsv files: 
```bash
cat file.tsv | head -1 | tr "\t" "\n" | wc -l  
csvcut -n -t  file.tsv (from csvkit)
awk '{print NF; exit}' file.tsv
awk -F "\t" 'NR == 1 {print NF}' file.tsv
```

#### combine info to the fasta header

[from biostar post](https://www.biostars.org/p/212379/#212393)
```bash
cat myfasta.txt 
>Blap_contig79
MSTDVDAKTRSKERASIAAFYVGRNIFVTGGTGFLGKVLIEKLLRSCPDVGEIFILMRPKAGLSI
>Bluc_contig23663
MSTNVDAKARSKERASIAAFYVGRNIFVTGGTGFLGKVLIEKLLRSCPDVGEIFILMRPKAGLSI
>Blap_contig7988
MSTDVDAKTRSKERASIAAFYVGRNIFVTGGTGFLGKVLIEKLLRSCPDVGEIFILMRPKAGLSI
>Bluc_contig1223663
MSTNVDAKARSKERASIAAFYVGRNIFVTGGTGFLGKVLIEKLLRSCPDVGEIFILMRPKAGLSI

cat my_info.txt 
info1
info2
info3
info4

paste <(cat my_info.txt) <(cat myfasta.txt| paste - - | cut -c2-) | awk '{printf(">%s_%s\n%s\n",$1,$2,$3);}'
>info1_Blap_contig79
MSTDVDAKTRSKERASIAAFYVGRNIFVTGGTGFLGKVLIEKLLRSCPDVGEIFILMRPKAGLSI
>info2_Bluc_contig23663
MSTNVDAKARSKERASIAAFYVGRNIFVTGGTGFLGKVLIEKLLRSCPDVGEIFILMRPKAGLSI
>info3_Blap_contig7988
MSTDVDAKTRSKERASIAAFYVGRNIFVTGGTGFLGKVLIEKLLRSCPDVGEIFILMRPKAGLSI
>info4_Bluc_contig1223663
MSTNVDAKARSKERASIAAFYVGRNIFVTGGTGFLGKVLIEKLLRSCPDVGEIFILMRPKAGLSI

```

#### count how many columns in a tsv file

```bash
cat file.tsv | head -1 | tr "\t" "\n" | wc -l  

##(from csvkit)
csvcut -n -t file.

## emulate csvcut -n -t
less files.tsv | head -1| tr "\t" "\n" | nl

awk -F "\t" 'NR == 1 {print NF}' file.tsv
awk '{print NF; exit}'
```
#### change fasta header

see https://www.biostars.org/p/53212/

The fasta header is like `>7 dna:chromosome chromosome:GRCh37:7:1:159138663:1`
convert to `>7`: 

```bash
cat Homo_sapiens_assembly19.fasta | gawk '/^>/ { b=gensub(" dna:.+", "", "g", $0); print b; next} {print}' > Homo_sapiens_assembly19_reheader.fasta
```
### mkdir and cd into that dir shortcut

```bash
mkdir blah && cd $_
```
### cut out columns based on column names in another file

http://crazyhottommy.blogspot.com/2016/10/cutting-out-500-columns-from-26g-file.html

```bash
#! /bin/bash

set -e
set -u
set -o pipefail

#### Author: Ming Tang (Tommy)
#### Date 09/29/2016
#### I got the idea from this stackOverflow post http://stackoverflow.com/questions/11098189/awk-extract-columns-from-file-based-on-header-selected-from-2nd-file

# show help
show_help(){
cat << EOF
  This is a wrapper extracting columns of a (big) dataframe based on a list of column names in another
  file. The column names must be one per line. The output will be stdout. For small files < 2G, one 
  can load it into R and do it easily, but when the file is big > 10G. R is quite cubersome. 
  Using unix commands on the other hand is better because files do not have to be loaded into memory at once.
  e.g. subset a 26G size file for 700 columns takes around 30 mins. Memory footage is very low ~4MB.

  usage: ${0##*/} -f < a dataframe  > -c < colNames> -d <delimiter of the file>
        -h display this help and exit.
		-f the file you want to extract columns from. must contain a header with column names.
		-c a file with the one column name per line.
		-d delimiter of the dataframe: , or \t. default is tab.  
		
		e.g. 
		
		for tsv file:
			${0##*/} -f mydata.tsv -c colnames.txt -d $'\t' or simply ommit the -d, default is tab.
		
		for csv file: Note you have to specify -d , if your file is csv, otherwise all columns will be cut out.
			${0##*/} -f mydata.csv -c colnames.txt -d ,
        
EOF
}

## if there are no arguments provided, show help
if [[ $# == 0 ]]; then show_help; exit 1; fi

while getopts ":hf:c:d:" opt; do
  case "$opt" in
    h) show_help;exit 0;;
    f) File2extract=$OPTARG;;
    c) colNames=$OPTARG;;
    d) delim=$OPTARG;;
    '?') echo "Invalid option $OPTARG"; show_help >&2; exit 1;;
  esac
done
	

## set up the default delimiter to be tab, Note the way I specify tab 

delim=${delim:-$'\t'}

## get the number of columns in the data frame that match the column names in the colNames file.
## change the output to 2,5,6,22,... and get rid of the last comma  so cut -f can be used
 
cols=$(head -1 "${File2extract}" | tr "${delim}" "\n" | grep -nf "${colNames}" | sed 's/:.*$//' | tr "\n" "," | sed 's/,$//')

## cut out the columns 
cut -d"${delim}" -f"${cols}" "${File2extract}"
```
or use [csvtk](https://github.com/shenwei356/csvtk) from Shen Wei:  

```bash
csvtk cut -t -f $(paste -s -d , list.txt) data.tsv
```
#### merge all bed files and add a column for the filename.

```bash
awk '{print $0 "\t" FILENAME}' *bed 
```

### add or remove chr from the start of each line

```bash
# add chr
sed 's/^/chr/' my.bed

# or
awk 'BEGIN {OFS = "\t"} {$1="chr"$1; print}'

# remove chr
sed 's/^chr//' my.bed
```
### check if a tsv files have the same number of columns for all rows

```bash
awk '{print NF}' test.tsv | sort -nu | head -n 1
```

### Parallelized samtools mpileup 

https://www.biostars.org/p/134331/

```bash
BAM="yourFile.bam"
REF="reference.fasta"
samtools view -H $BAM | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d 100000 -uf $REF -r \"{}\" $BAM | bcftools call -cv > \"{}\".vcf"
```
### convert multiple lines to a single line

This is better than `tr "\n" "\t"` because somtimes I do not want to convert the last newline to tab.

```bash
cat myfile.txt | paste -s 
```

### merge multiple files with same header by keeping the header of the first file
I usually do it in R, but like the quick solution.

https://stackoverflow.com/questions/16890582/unixmerge-multiple-csv-files-with-same-header-by-keeping-the-header-of-the-firs

```bash
awk 'FNR==1 && NR!=1{next;}{print}' *.csv 

# or

awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' file*.txt >all.txt
```

### insert a field into the first line

```bash
cut -f1-4 F5.hg38.enhancers.expression.usage.matrix | head
CNhs11844	CNhs11251	CNhs11282	CNhs10746
chr10:100006233-100006603	1	0	0
chr10:100008181-100008444	0	0	0
chr10:100014348-100014634	0	0	0
chr10:100020065-100020562	0	0	0
chr10:100043485-100043744	0	0	0
chr10:100114218-100114567	0	0	0
chr10:100148595-100148922	0	0	0
chr10:100182422-100182522	0	0	0
chr10:100184498-100184704	0	0	0

sed '1 s/^/enhancer\t/' F5.hg38.enhancers.expression.usage.matrix | cut -f1-4 | head
enhancer	CNhs11844	CNhs11251	CNhs11282
chr10:100006233-100006603	1	0	0
chr10:100008181-100008444	0	0	0
chr10:100014348-100014634	0	0	0
chr10:100020065-100020562	0	0	0
chr10:100043485-100043744	0	0	0
chr10:100114218-100114567	0	0	0
chr10:100148595-100148922	0	0	0
chr10:100182422-100182522	0	0	0
chr10:100184498-100184704	0	0	0

```
### extract PASS calls from vcf file

```
cat my.vcf | awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' > my_PASS.vcf

```

### replace a pattern in a specific column

```
## column5 
awk '{gsub(pattern,replace,$5)}1' in.file

## http://bioinf.shenwei.me/csvtk/usage/#replace
csvtk replace -f 5 -p pattern -r replacement 

```
### move a process to a screen session

https://www.linkedin.com/pulse/move-running-process-screen-bruce-werdschinski/

```
1. Suspend: Ctrl+z
2. Resume: bg
3. Disown: disown %1
4. Launch screen
5. Find pid: prep BLAH
6. Reparent process: reptyr ###
```

### count uinque values in a column and put in a new 

https://www.unix.com/unix-for-beginners-questions-and-answers/270526-awk-count-unique-element-array.html

```
# input
blabla_1 A,B,C,C
blabla_2 A,E,G
blabla_3 R,Q,A,B,C,R,Q

# output
blabla_1 3
blabla_2 3
blabla_3 5


awk '{split(x,C); n=split($2,F,/,/); for(i in F) if(C[F[i]]++) n--; print $1, n}' file

```

### get the promoter regions from a gtf file

https://twitter.com/David_McGaughey/status/1106371758142173185

Create TSS bed from GTF in one line: 
```bash
zcat gencode.v29lift37.annotation.gtf.gz | awk '$3=="gene" {print $0}' | grep protein_coding | awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+1} else {print $1, $5-1, $5}}' > tss.bed
```
or 5kb flanking tss

```bash
zcat gencode.v29lift37.annotation.gtf.gz | awk '$3=="gene" {print $0}' | grep protein_coding | awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+5000} else {print $1, $5-5000, $5}}' > promoters.bed
```
caveat: some genes are at the end of the chromosomes, add or minus 5000 may go beyond the point, use [`bedtools slop`](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html) with a genome size file to avoid that.

download `fetchChromSizes` from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

```bash
fetchChromSizes hg19 > chrom_size.txt

zcat gencode.v29lift37.annotation.gtf.gz | awk '$3=="gene" {print $0}' |  awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+1} else {print $1, $5-1, $5}}' | bedtools slop -i - -g chrom_size.txt -b 5000 > promoter_5kb.bed
```

### reverse one column of a txt file

reverse column 3 and put it to column5
```bash
awk -v OFS="\t" '{"echo "$3 "| rev" | getline $5}{print $0}' 

#or use perl reverse second column
perl -lane 'BEGIN{$,="\t"}{$rev=reverse $F[2];print $F[0],$F[1],$rev,$F[3]}
```
