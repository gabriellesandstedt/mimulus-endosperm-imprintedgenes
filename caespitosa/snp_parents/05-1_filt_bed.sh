awk 'BEGIN{OFS="\t"}
$11 ~ /^[01][\/|][01]/ && $12 ~ /^[01][\/|][01]/ {
    split($11,a,":"); gtA=a[1]
    split($12,b,":"); gtB=b[1]

    if (gtA=="0/0" || gtA=="0|0") alleleA=$6
    else if (gtA=="1/1" || gtA=="1|1") alleleA=$7
    else next

    if (gtB=="0/0" || gtB=="0|0") alleleB=$6
    else if (gtB=="1/1" || gtB=="1|1") alleleB=$7
    else next

    print $1,$2,$3,alleleA">"alleleB
}' til.bed > diagnostic_snps_tilsp_SOP_LVR_caesref.bed


#### sop12 is sp A, LVR sp B
### utc1 is sp A, twn is spB
