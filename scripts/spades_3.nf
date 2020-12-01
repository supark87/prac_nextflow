//params.reads='/Users/subinpark/fingerprinting/D0/*_{R1,R2}_001.fastq.gz'
params.reads='/Users/subinpark/NeST_CDC/NeST/input2/*_{R1,R2}_001.fastq.gz'
params.references='/Users/subinpark/fingerprinting/reference_genes/cpmp_NC_004325.fasta'
Channel
    .fromFilePairs(params.reads)
    .into {reads_ch; reads_ch2}

//Channel
//   .fromPath(params.cpmp)
//params.alignments="/Users/subinpark/prac_nextflow/result_tset_spades_2/muscle_output/*.fasta"
//params.alignments='/Users/subinpark/fingerprinting/psf_cpmp_no_prediction/blastoutput/seqrecords/pergene_seqrecords/muscle_output'

// Channel
//      .fromPath(params.alignments)
//      .into {alignments_ch; call_genotype_ch}
params.outdir="individual_whole_test"

//type, variable, channel for input and output

process fastqc {
    container 'my-image-spades'
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_out"
    
    input:
    tuple val(sample_id), file(reads_file) from reads_ch 

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}
process bbduk{
    container 'my-image-spades'

    tag "$sample_id"
    publishDir "${params.outdir}/bbduk_out", mode:'copy'

    input:
    tuple val(sample_id), file(reads_file) from reads_ch2

    output:
    tuple val(sample_id), file("${sample_id}_clean_{1,2}.fastq") into skesa_input
    file "${sample_id}.stats.txt"

    script:
    """
    bbduk.sh -Xmx1g in=${reads_file[0]} in2=${reads_file[1]} out=${sample_id}_clean_1.fastq out2=${sample_id}_clean_2.fastq qtrim=rl trimq=30 minlength=50\
    stats=${sample_id}.stats.txt
    
    """

}


process spades{
    container 'pegi3s/spades'

    tag "$sampe_id"
    publishDir "${params.outdir}/spades_out", mode:'copy'
    //errorStrategy 'ignore'
    errorStrategy 'retry'

    input:
    tuple val(sample_id),file(reads_file) from skesa_input

    output:
    tuple val(sample_id), file("${sample_id}_scaffolds.fasta") into spade_output,spade_output_for_prod
     
    script:
    """
    spades.py -k 21,33,55,77,99 --careful --pe1-1 ${reads_file[0]} --pe1-2 ${reads_file[1]}  -o ${sample_id} 
    cp ${sample_id}/scaffolds.fasta ${sample_id}_scaffolds.fasta
    """
}

// process run_prodigal{
//     tag "$sample_id"
//     publishDir "${params.outdir}/prod_fasta",mode:'copy'

//     input:
//     set val(sample_id), file(assembly) from spade_output_for_prod

//     output:
//     set val(sample_id), file("${assembly}_prod.fasta") into prod_fasta

//     when:
//     assembly.size() > 0

//     script:
//     """
//     prodigal -i ${assembly} -p meta -d ${assembly}_prod.fasta
//     """
// }

process makeblastdb{
    container 'my-image-spades'

    tag "$sample_id"
    publishDir "${params.outdir}/blast_db",mode:'copy'

    input:
    //set val(sample_id),file(assembly) from prod_fasta
    set val(sample_id), file(assembly) from spade_output_for_prod

    output:
    val "${sample_id}_db" into blastdb_name
    path "${sample_id}_dir" into blastdb_dir

    when:
    assembly.size() > 0
    
    script:
    """
    mkdir ${sample_id}_dir
    makeblastdb -in ${assembly} -dbtype nucl -out ${sample_id}_dir/${sample_id}_db
    """


}

process blastrun{
    container 'my-image-spades'

    tag "$sample_id"
    publishDir "${params.outdir}/blast_out_files",mode:'copy'
    
    input:
    path cpmp from params.references
    val sample_id_db from blastdb_name 
    path dbdir from blastdb_dir

    output:
    file ("${sample_id_db}.blast") into blastout_ch
    file("${sample_id_db}.fasta") into sequences
    val "${sample_id_db}" into seq_name
    file("${sample_id_db}.fasta") into cpmp_all
    
    script:


    """
    blastn -query $cpmp -db $dbdir/$sample_id_db  -outfmt "6 qseqid sseqid pident qlen qstart qend sstart send sseq" -out ${sample_id_db}.blast
    cat ${sample_id_db}.blast | awk '\$8-\$7 > 300 {print ">" \$2 "_${sample_id_db}_" "\\n" \$9}' > ${sample_id_db}.fasta
    
    """
}

//params.allcpmp=sequences.collectFile()
      
process alginment{
    container 'my-image-spades'

    tag "$sample_id"
    publishDir "${params.outdir}/muscle_output",mode:'copy'

    input:
    file('*') from sequences.collect()
    val sample_id_db from seq_name

    output:
    file("cpmp.align.fasta") into (alignments_ch, call_genotype_ch)

    script:
    """
    cat * > cpmp_all.fasta
    muscle -in cpmp_all.fasta -out cpmp.align.fasta
    """

 }


process unique_allele{
    container='biopython/biopython'
    publishDir "${params.outdir}/allele_subtype",mode:'copy'
    input:
    file(alignments) from alignments_ch
    output:
    file("${alignments}_uniq.fasta") into (unique_allele_ch1,unique_allele_ch2,unique_allele_ch3)
    file("${alignments}_uniq.fasta") into unique_allele_ch4
    stdout into corename_ch

    script:
    """
    #!/usr/bin/python3

    from Bio import AlignIO
    import os
    import sys

    corename="$alignments".split(".")[0]
    print(corename.rstrip())
    count=1
    alignment=AlignIO.read("$alignments","fasta")
    dic1={}
    for record in alignment :
        dic1[str(record.seq)]=1
    for key1 in dic1:
        with open("$alignments"+"_"+"uniq.fasta",'a') as ofile:
            ofile.write(">"+corename+"_"+str(count)+"\\n"+key1+"\\n")
        count=count+1

    """

 }

process assign_allele{
    container='biopython/biopython'
    publishDir "${params.outdir}/Genomeprofile",mode:'copy'

    input:
    file(allele_subtype) from unique_allele_ch2
    file(alignment) from call_genotype_ch

    output:
    file("Genome_profile.csv") into alignment_match_genotype_ch
    
    script:
"""
#!/usr/bin/python3
from Bio import AlignIO
import os
import sys
dic1={}
dic2={}
dic3={}
def Parse(filename,seqs):
    file = open(filename)
    seqs={}
    name = ''
    for line in file:
        line = line.rstrip()  
        if line.startswith('>'):     
            name=line.replace('>',"") 
            seqs[name] = ''          
        else:
            seqs[name] = seqs[name] + line
    return seqs   
dic1=dict(dic1,**Parse("$allele_subtype",dic1))

corename="$alignment".split("_")[0]
alignment=AlignIO.read("$alignment","fasta")

for record in alignment :    
    #record.id1=record.id.split("_")[::-1][0]
    dic2.setdefault(str(record.id),[]).append(str(record.seq))
for key1 in dic1:
    for key2 in dic2:
        if dic1.get(key1) in dic2.get(key2):
            dic3.setdefault(str(key2),[]).append(str(key1))

for key3 in dic3:
    with open('Genome_profile.csv','a') as ofile:
        ofile.write(str(key3)+"\\t"+str(dic3.get(key3)[0])+"\\n")
        ofile.close()
"""

}

process macse_alignment{
    container 'ranwez/omm_macse:v10.02'
    publishDir "${params.outdir}/macse_output",mode:'copy'

     input: 
     file(unique_allele_subtype) from unique_allele_ch4
     val(core_name) from corename_ch

     output:
     //file('*') into macse_outfiles
     path ('*') into macse_path

     script:
     """
     /OMM_MACSE/S_OMM_MACSE_V10.02.sh --out_dir macse_out  --in_seq_file "$unique_allele_subtype"  --alignAA_soft MAFFT  --out_file_prefix ${core_name}
     """

 }

process call_genotype{
    container='thibautjombart/adegenet_testing:latest'
    publishDir "${params.outdir}/genotype",mode:'copy'
    input:
    file(alignments) from unique_allele_ch1
    output:
    file("cpmp_genotype_profile.csv") into genotype_profile_ch
    file("cpmp_snp_density.png")
    script:
    """
    #!/usr/bin/Rscript --vanilla
    library(ape)
    library(adegenet)
    library(poppr)
    cpmp<-fasta2DNAbin("$alignments")
    cpmp_snp<-as.matrix(cpmp)
    png("cpmp_snp_density.png")
    snpposi.plot(cpmp_snp,codon=TRUE)
    dev.off()
    obj_cpmp = DNAbin2genind(cpmp_snp, polyThres=0.1)
    mlg(obj_cpmp)
    genind2genalex(obj_cpmp, filename = "cpmp_genotype_profile.csv", quiet = FALSE,
    pop = NULL, allstrata = TRUE, geo = FALSE, geodf = "xy",overwrite = TRUE,
    sep = ",") 
    df_cpmp=genind2df(obj_cpmp)

    """
}

process alignment_match_genotype {
    container='biopython/biopython'
    publishDir "${params.outdir}/allele_assigne_new_name",mode:'copy'
    input:
    file(genotype_profile) from alignment_match_genotype_ch
    file(allele_profile) from genotype_profile_ch
    output:
    file("new_name_allele_call.csv")
    file("new_name_allele.csv")
    script:

 """
#!/usr/bin/python3
import os
import pandas as pd
import numpy as np

def f(row):
    if row['D0_new_name'] == row['DF_new_name']:
        val = 0
    else:
        val = -1
    return val

D0_DF=pd.read_csv("$genotype_profile",header=None,sep="\t")
D0_DF['sample']=D0_DF[0].apply(lambda x:x.split("_")[::-1][4])
D0=D0_DF[D0_DF['sample'].apply(lambda x:x[6:8])=="00"]
DF=D0_DF[D0_DF['sample'].apply(lambda x:x[6:8]!="00")]
D0['sample_id']=D0['sample'].apply(lambda x:x[4:6])+"_"+D0['sample'].apply(lambda x:x[8:13])
DF['sample_id']=DF['sample'].apply(lambda x:x[4:6])+"_"+DF['sample'].apply(lambda x:x[8:13])
D0.columns=['all_name','D0_allele','sample','sample_id']
DF.columns=['all_name','DF_allele','sample','sample_id']

combine=pd.merge(D0,DF,on='sample_id',how='outer')
combine=pd.merge(D0,DF,on='sample_id',how='outer')
combine.columns=['D0_sample_id','D0_allele','D0_sample','sample_id','DF_sample_id','DF_allele','DF_sample']
combine['DayF']=combine.DF_sample.str[6:8]
combine_sub=combine[['sample_id','D0_sample','D0_allele','DayF','DF_allele']]
df=combine_sub[combine_sub.DayF!="xx"].sort_values(by='sample_id')

allele_profile=pd.read_csv("$allele_profile")
allele_profile.columns=allele_profile.iloc[1]
allele_profile=allele_profile.drop(0,axis=0).drop(1,axis=0)
allele_profile=allele_profile.drop("Pop",axis=1)
allele_profile['id']= allele_profile[allele_profile.columns[1:]].apply(
    lambda x: ','.join(x.dropna().astype(str)),
    axis=1
)
allele_profile_table=pd.DataFrame(allele_profile.groupby('id')['Ind'].agg(lambda x:' '.join(x.unique())))
allele_profile_table=allele_profile_table.reset_index()
allele_profile_table=allele_profile_table.reset_index()
allele_profile_table.columns=['unique_id','id','Ind']
allele_profile_table['new_name']=allele_profile_table.apply(lambda row: "cpmp_" + str(row['unique_id']), axis=1)
allele_profile_table.to_csv("new_name_allele.csv")

name_allele_map=pd.merge(allele_profile,allele_profile_table,on='id')[['Ind_x','new_name']]
name_allele_map2=name_allele_map
name_allele_map.columns=['D0_allele','D0_new_name']
merge1=pd.merge(df,name_allele_map,on='D0_allele',how='left')
name_allele_map2.columns=['DF_allele','DF_new_name']
merge2=pd.merge(merge1,name_allele_map2,on='DF_allele',how='left')
merge2=merge2[['sample_id','D0_new_name','DayF','DF_new_name']]

merge2['recrudescence'] = merge2.apply(f, axis=1)
merge2.to_csv("new_name_allele_call.csv")

 """

 }

