1) Download data from NCBI-SRA

ERX535560 adult      ERR577157 ad2     
ERX535563 adult      ERR577160 ad3   
ERX535561 egg       ERR577158	egg2   
ERX535360 egg       ERR576954    egg1   
ERX535363 juvenile    ERR576957    juv1   
ERX535362 juvenile     ERR576956    juv1   
ERX535371 metacercaria   ERR576965   met0hr_3   
ERX535370 metacerc.   ERR576964    met0hr_2   


```bash
prefetch ERR*
```
```bash
for i in ERR*; do  
fastq-dump -I --split-3 $i
done
``` 

2) Quality control

```bash
for i in ERR*.fastq; do  
fastqc --noextract $i
done
```

3) Trimming

```bash
for i in ERR*.fastq; do    
java -jar ./Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 $i -baseout ./trim/$i ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:7 LEADING:3 TRAILING:3 MAXINFO:25:0.4 MINLEN:25
done
```

4) Quality control


```bash
for i in ERR*.fastq; do  
fastqc --noextract $i
done
```

5) DE analysis
  5.1 kallisto+sleuth (https://www.ncbi.nlm.nih.gov/pubmed/28581496)
 
kallisto 
```bash  
kallisto quant -i repeats -o ERR* -t 8 -b 100 ./trim/ERR*_1P.fastq ./trim/ERR*_2P.fastq
```

sleuth

info.txt:
sample	condition
576954	egg
576956	juv
576957	juv
576964	metacerc
576965	metacerc
577157	adult
577158	egg
577160	adult

```r
library(sleuth)
base_dir <- "~/kallisto_b"
sample_id <- dir(file.path(base_dir))
s2c <- read.table(file.path(base_dir, "info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample, condition)                 
kal_dirs <- file.path(base_dir, sample_id)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

so <- sleuth_prep(s2c, extra_bootstrap_summary = T)

so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

sleuth_live(so)
```
