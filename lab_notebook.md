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
results:  
see [kallisto_b](https://github.com/NickPanyushev/IB_Fasciola/tree/master/kallisto_b)

sleuth

info.txt:    
see [info.txt](https://github.com/NickPanyushev/IB_Fasciola/blob/master/kallisto_b/info.txt)  

code:
see [sleuth.R](https://github.com/NickPanyushev/IB_Fasciola/blob/master/sleuth.R)
results:
see [sleuth_res](https://github.com/NickPanyushev/IB_Fasciola/tree/master/sleuth_res)

  5.1 TEtools (TEcount+DeSeq2)
 rosette file:  
 see [same.csv](https://github.com/NickPanyushev/IB_Fasciola/blob/master/same.csv)  
 
 TEcount  
 ```bash
  python3 -u TEcount.py -rosette same.csv -column 3 -TE_fasta f_hepatica-families.fa -count all_same_3 -RNA /Johnny/skalon/fasciola/trim/ERR577157_filtered5_1P.fastq /Johnny/skalon/fasciola/trim/ERR577160_filtered_1P.fastq /Johnny/skalon/fasciola/trim/ERR577158_filtered_1P.fastq /Johnny/skalon/fasciola/trim/ERR576954_filtered_1P.fastq /Johnny/skalon/fasciola/trim/ERR576957_filtered_1P.fastq /Johnny/skalon/fasciola/trim/ERR576956_filtered_1P.fastq /Johnny/skalon/fasciola/trim/ERR576965_filtered_1P.fastq /Johnny/skalon/fasciola/trim/ERR576964_filtered_1P.fastq -RNApair /Johnny/skalon/fasciola/trim/ERR577157_filtered5_2P.fastq /Johnny/skalon/fasciola/trim/ERR577160_filtered_2P.fastq /Johnny/skalon/fasciola/trim/ERR577158_filtered_2P.fastq /Johnny/skalon/fasciola/trim/ERR576954_filtered_2P.fastq /Johnny/skalon/fasciola/trim/ERR576957_filtered_2P.fastq /Johnny/skalon/fasciola/trim/ERR576956_filtered_2P.fastq /Johnny/skalon/fasciola/trim/ERR576965_filtered_2P.fastq /Johnny/skalon/fasciola/trim/ERR576964_filtered_2P.fastq -bowtie2 -insert 100 | tee all_same_3.log
  ```
  results:  
  see [all_same_3.csv](https://github.com/NickPanyushev/IB_Fasciola/blob/master/all_same_3.csv)  
  
  DeSeq2.R 
  ```bash
   Rscript DeSeq2.R --args --FDR_level=0.05 --count_column=3 --count_file=\"all_same_3.csv\" experiment_formula=\"sample:replicant\" --sample_names=\"adult:1,adult:2,egg:1,egg:2,juv:1,juv:2,met:1,met:2\" --outdir=\"DeSeq2_res\"
```
code:    
see [DeSeq2.R](https://github.com/NickPanyushev/IB_Fasciola/blob/master/DeSeq2.R)  
results:   
see [DeSeq2_res](https://github.com/NickPanyushev/IB_Fasciola/tree/master/DeSeq2_res)
