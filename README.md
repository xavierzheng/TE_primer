# Objective

This tutorial outlines a potential method for designing primer sets to specifically amplify a transposon.



# Methods

## Step 0: Create a Unique Folder for Each TE

Before starting any analysis, I recommend creating a unique folder for each TE, named after its unique ID.



## Step 1: Generate Multiple Primer Sets with Primer3

1. Go to primer3: [Primer3 Input](https://primer3.ut.ee/)

2. Copy and paste the transposon sequence from a genome browser, for example, [BoA Jbrowse](http://bcg.abrc.sinica.edu.tw/JBrowse/BoA/)
3. Modify these parameters:
   * **Product size:** 150-200
   * **Max 3' stability:** 12
   * **Number To Return:** 100
   * **Primer size:** 18-20-21
   * **Tm:** 59-60-61 (different: 0.5)
   * **Max self complement (or PRIMER_MAX_SELF_ANY):** 3~5, lower is better
   * **Max 3' complement (or PRIMER_MAX_SELF_END):** 2
   * **Poly X (or PRIMER_MAX_POLY_X):** 3
   * **GC content:** 30-70 or 40-60
   * **Salt Correction Formula:** SantaLucia1998 
   * **------ below parameters are depent on the qPCR kit-------**
   * **Concentration of divalent cations:** 2
   * **Concentration of dNTPs:** 0.8
   * **Annealing Oligo Concentration:** 200

4. Copy all primer results in the following format into a single file:
   * Primer set name
   * Underline
   * F/R (forward/reverse)

5. For example:

```tex
>LTR1_F
CGGAGAATGATGGTGGTCTA
>LTR1_R
ACTTGTTACGCTCCCTTGTA
>LTR2_F
AGAAGAGTCAGCTACACCAC
>LTR2_R
GGCTCTCAATCCGCAAATAG
>LTR3_F
CACCAACTTCACCAACAACA
>LTR3_R
TGCACTTCTCCATGACTTCT
```



## Step 2: Check Primer Specificity with Blastn

1. Before running the code, ensure you have installed `blast 2.13.0`
2. Run `blastn` locally using `transform_PrimerBlastN.sh`

```bash
transform_PrimerBlastN.sh <primer file> <output result file name>
```



## Step 3: Select the Primer Set with the Highest Specificity

At this stage, we'll select the primer set with 100% coverage, <2 mismatches, and 0 gaps from the step 2 `transform_PrimerBlastN.sh` output. Then, we'll select amplicons less than 500 bp in size as larger fragments are more challenging to amplify under qPCR conditions.

1. Before using this R script, make sure you have installed `R (version >= 3.5.0)` and the `tidyverse` library
2. Run the following scrip:

```bash
Rscript transform_PrimerBlastTobed.R <output result file name from step 2>
```

3. This script will generate 3 files, including `primer_target_site.bed`, `potential_amplicon.bed`,  and `potential_amplicon_summary.txt`

   * `primer_target_site.bed`: After selecting the blastn result with 100% coverage with <2 mismatch and 0 gap, this script will generate the selected results in bed 6 column format

     * For example: the 5th column is the number of mismatch base

       ```text
       BoA_C1  7251898 7251918 LTR12_F 1       +
       BoA_C1  7251917 7251939 LTR91_R 1       -
       BoA_C1  7262066 7262086 LTR53_R 0       -
       BoA_C1  7262068 7262088 LTR51_F 0       +
       BoA_C1  7262289 7262311 LTR73_F 1       +
       BoA_C1  7262354 7262375 LTR84_R 0       -
       BoA_C1  7262355 7262376 LTR57_F 0       +
       ```

       

   * `potential_amplicon.bed`: In this file, only the amplicons less than 500bp are recorded in bed 6 column format

     * For example: the 5th column is the amplicon size

       ```text
       BoA_C1  49785647        49785826        LTR76   179     +
       BoA_C2  29261759        29261938        LTR76   179     +
       BoA_C2  31447061        31447240        LTR76   179     +
       BoA_C2  57624045        57624224        LTR76   179     +
       BoA_C3  1751492 1751671 LTR76   179     +
       BoA_C3  34230729        34230908        LTR76   179     +
       BoA_C3  35347267        35347446        LTR76   179     +
       ```

   * `potential_amplicon_summary.txt`: In this file, the number of potential amplicons in specific genome are shown.

     * For example. "GRP" column indicates the primer set name. "COUNT" column indicates the number of potential amplicons. in this example, LTR12 and LTR14 have only one amplicon. This is the candidate primer set which can specifically amplify the TE

       ```text
       GRP     COUNT
       LTR1    31
       LTR10   75
       LTR100  29
       LTR11   86
       LTR12   1
       LTR13   71
       LTR14   1
       LTR15   2
       LTR16   5
       ```

       

## Step 4: Identify the Primer Set That Can Uniquely Amplify the TE

Now we'll parse the results from Step 3 and generate three files (`specific_set.txt`, `specific_primer_target_site.bed`, and `specific_potential_amplicon.bed`) using this script:

```bash
bash transform_PrimerBlastToSpecific.sh
```

* `specific_set.txt`: this file record the primer set which have only one amplicon within specific genome. Then, I add one underline in primer set name 

  * For example:

    ```tex
    LTR12_
    LTR14_
    LTR34_
    LTR46_
    LTR5_
    LTR50_
    ```

* `specific_primer_target_site.bed`: This file records the primer target sites. It is different from the results of `primer_target_site.bed` in step 3, because in this file, we only record the primer sets listed in `specific_set.txt`

  * For example:

    ```text
    BoA_C1  460971  460992  LTR1572lot68_F  1       -
    BoA_C1  2270197 2270218 LTR1572lot68_F  1       -
    BoA_C1  7249426 7249447 LTR1572lot86_R  0       -
    BoA_C1  7250521 7250541 LTR1572lot5_R   1       -
    BoA_C1  7250704 7250724 LTR1572lot34_R  0       -
    BoA_C1  7251462 7251482 LTR1572lot46_F  1       +
    BoA_C1  7251596 7251616 LTR1572lot14_F  0       +
    ```

* `specific_potential_amplicon.bed`: This file records the primer amplicon. It is different from the results of `potential_amplicon.bed` in step 3, because in this file, we only record the primer amplicon listed in `specific_set.txt`

  * For example:

    ```text
    BoA_C3  72536631        72536792        LTR1572lot68    161     +
    BoA_C3  72538071        72538228        LTR1572lot5     157     +
    BoA_C3  72538214        72538413        LTR1572lot34    199     +
    BoA_C3  72539287        72539441        LTR1572lot14    154     +
    BoA_C3  72539141        72539311        LTR1572lot95    170     +
    BoA_C3  72539328        72539520        LTR1572lot56    192     +
    ```



## Step 5: Add TE ID to Primer Set Name

The results from Step 4 might not include the TE ID in the primer set name. Thus, in this step, we add the TE ID to the primer set name using the following script:

```bash
bash transform_PrimerBlastAddPrefix.sh <TE ID>
```

Ensure your folder name matches your TE ID. After this, manually check the results. If everything looks fine, remember to remove the temp files:

```bash
rm <TE ID>/specific_potential_amplicon.bed.temp
rm <TE ID>/specific_primer_target_site.bed.temp
```



## Step 6: Visualize Your Results

Import the results of `specific_potential_amplicon.bed` and `specific_primer_target_site.bed` from Step 4 into a genomics visualization tool such as IGV, Jbrowse, or CLC Genomic Workbench. In combination with a WGS bam file, you can visually assess your primer sets for any potential issues.
