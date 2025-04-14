# IBD SMR # 
pattern="(.+)_gwas\.ma"
cat GTEx.txt | while read tissue
do
    ls GWASIBD | while read id
    do 
            if [[ $id =~ $pattern ]]; then
                outid=${BASH_REMATCH[1]}
                echo $outid
            else
                echo "未匹配到模式"
            fi
        ./smr-1.3.1-linux-x86_64/smr-1.3.1 --bfile ./EUR --gwas-summary ./GWASIBD/"$id" --beqtl-summary ./"$tissue"/"$tissue" --out ./"$tissue"_ibdout/"$outid" --thread-num 10
    done    
done



