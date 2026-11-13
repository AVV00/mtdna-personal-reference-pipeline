#!/bin/bash
start=$(date +%s) #脚本启始时间
echo "开始时间："`date` 

#设置文件的路�?
sample_dir="/mnt/data2/lishengjing/2_o_10_HCC/DNA/"
result_dir="/mnt/data2/lishengjing/2_o_10_HCC/DNA/"


cd $result_dir
for i in $(cat ${sample_dir}/sample_name.txt)
do
    #mkdir ${i}
    echo "
    #!/bin/bash
    cd $result_dir/${i}
    mkdir tmp
    Rscript /mnt/hdd/usr/xiefanfan/work/NUMT/before/script/modified.mt.genome.R
    cat /mnt/hdd/ref_data/genomes/human/hg38/hg38_nochrM.fasta  modified.${i}.fa > ${i}.hg38.mt.fa
    bwa index ${result_dir}/${i}/${i}.hg38.mt.fa
    bwa mem -t 10 -M -R '@RG\tID:${i}\tLB:mtDNA\tPL:ILLUMINA\tSM:${i}' ${sample_dir}/${i}/${i}.hg38.mt.fa ${sample_dir}/${i}/${i}_R1.good.fq.gz ${sample_dir}/${i}/${i}_R2.good.fq.gz 2>>${i}_analysis.log |samtools view -bS - -o ${i}.raw.bam && 
    java -Xmx10g -Djava.io.tmpdir=`pwd`/tmp -jar /mnt/hdd/softs/picard/picard-tools-1.81/SortSam.jar I=${i}.raw.bam O=${i}.sort.picard.bam SO=coordinate VALIDATION_STRINGENCY=SILENT 2>>${i}_analysis.log   && 
    java -Xmx10g -Djava.io.tmpdir=`pwd`/tmp -jar /mnt/hdd/softs/picard/picard-tools-1.81/BuildBamIndex.jar INPUT=${i}.sort.picard.bam OUTPUT=${i}.sort.picard.bai VALIDATION_STRINGENCY=SILENT 2>>${i}_analysis.log && 
    java -Xmx10g -Djava.io.tmpdir=`pwd`/tmp -jar /mnt/hdd/softs/picard/picard-tools-1.81/MarkDuplicates.jar REMOVE_DUPLICATES=true I=${i}.sort.picard.bam O=${i}.rmdup.bam M=${i}.dup.metric VALIDATION_STRINGENCY=SILENT 2>>${i}_analysis.log  && 
    java -Xmx10g -Djava.io.tmpdir=`pwd`/tmp -jar /mnt/hdd/softs/picard/picard-tools-1.81/BuildBamIndex.jar INPUT=${i}.rmdup.bam  OUTPUT=${i}.rmdup.bai VALIDATION_STRINGENCY=SILENT 2>>${i}_analysis.log &&  
    perl /mnt/hdd/softs/mt_analysis/fetch_mt.pl ${i}.rmdup.bam |samtools view -bS - -o ${i}.mt.bam  2>>${i}_analysis.log && 
    samtools index ${i}.mt.bam && 
    samtools depth -m 0 ${i}.mt.bam > ${i}.mt.depth && 
    samtools view ${i}.mt.bam > ${i}.mt.sam  && 
    python2 /mnt/hdd/softs/DREEP/sam_to_pileup.py ${i}.mt.sam modified.${i}.fa > ${i}.old.pileup  && 
    perl /mnt/hdd/softs/DREEP/pileup_appended.pl ${i}.old.pileup > ${i}.pileup && 
    /usr/bin/perl /mnt/hdd/usr/lishengjing/Code/filter_and_summary_v3.pl -i ${i}.pileup -d 3 -q 30 -m 20 -r 143 -s 5  > ${i}.new.ssp 2>>${i}_analysis.log && 
    rm ${i}.old.pileup && 
    rm ${i}.mt.sam && 
    rm ${i}.pileup &&  
    rm -rf tmp
        " > ${i}/${i}_analysis.sh
done


#for循环并行运算
tmp_file="/tmp/$$.fifo"
mkfifo $tmp_file
exec 6<>$tmp_file   
rm -f $tmp_file
thread=12   #设置并行运算的数�?           
for ((i=0;i<thread;i++))    
do
    echo >&6
done



#第二步：运行脚本
#sh $code_file_dir/run_task.sh  
#!/bin/bash
function feth_file(){
    for file in $(cat ${sample_dir}/sample_name.txt)
    do  
    read -u6
   {
        df=$1"/"$file
        if [ -d $df ]
        then 
            cd $df
            bash *.sh
            cd ..
        fi
        echo >&6
   }&
    done
    wait        # 等待所有的后台子进程结�?    exec 6>&-   # 关闭fd6          
}
feth_file $result_dir  
rm -rf $result_dir/tmp


end=$(date +%s) #脚本终止时间
take=$(( end - start ))
h=$(( take / 3600 ))
m=$((( take - 3600 * h ) / 60 ))
s=$((( take - 3600 * h - 60 * m )))
echo "结束时间�?`date` 
echo "脚本运行完毕，用�?${h} 小时 ${m} 分钟 ${s} �?
