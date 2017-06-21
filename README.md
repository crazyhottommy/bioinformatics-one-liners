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

####bam2wig

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

####split large file by id/label/column

```bash
awk '{print >> $1; close($1)}' input_file
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
