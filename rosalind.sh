#!/bin/bash
# rosalind.sh

# dna0() {
#     # Fails edge cases
#     echo "$1" | grep -o . | sort | uniq -c | awk '{print $1}' ORS=' '
# }
# 
# dna1() {
#     local a=0
#     local c=0
#     local g=0
#     local t=0
#     for ((i = 0; i < ${#1}; i++)); do
#         [[ "${1:$i:1}" == 'A' ]] && ((++a)) && continue
#         [[ "${1:$i:1}" == 'C' ]] && ((++c)) && continue
#         [[ "${1:$i:1}" == 'G' ]] && ((++g)) && continue
#         [[ "${1:$i:1}" == 'T' ]] && ((++t)) && continue
#     done
#     echo "$a $c $g $t"
# }

dna() {
    local arr
    declare -A arr
    arr["A"]=0
    arr["C"]=0
    arr["G"]=0
    arr["T"]=0
    for ((i = 0; i < ${#1}; i++)); do
        ((++arr["${1:$i:1}"]))
    done
    echo "${arr["A"]} ${arr["C"]} ${arr["G"]} ${arr["T"]}"
}

test_dna() {
    local actual
    actual=$(dna "$1")
    if [ "$actual" != "$2" ]; then
        echo "failed dna=$1 expected=$2 actual=${actual}"
        exit 1
    fi
}

test_all_dna() {
    test_dna '' '0 0 0 0'
    test_dna 'A' '1 0 0 0'
    test_dna 'C' '0 1 0 0'
    test_dna 'G' '0 0 1 0'
    test_dna 'T' '0 0 0 1'
    test_dna 'ACGT' '1 1 1 1'
    test_dna 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC' '20 12 17 21'
    echo "DNA : PASS"
}

# test_all_dna



rna() {
    tr T U <<< "$1"
}

test_rna() {
    local actual
    actual=$(rna "$1")
    if [ "$actual" != "$2" ]; then
        echo "failed rna=$1 expected=$2 actual=${actual}"
        exit 1
    fi
}

test_all_rna() {
    test_rna '' ''
    test_rna 'A' 'A'
    test_rna 'C' 'C'
    test_rna 'G' 'G'
    test_rna 'T' 'U'
    test_rna 'GATGGAACTTGACTACGTAAATT' 'GAUGGAACUUGACUACGUAAAUU'
    echo "RNA : PASS"
}

# test_all_rna



revc() {
    tr ACGTacgt TGCAtgca <<< "$1" | rev
}

test_revc() {
    local actual
    actual=$(revc "$1")
    if [ "$actual" != "$2" ]; then
        echo "failed revc=$1 expected=$2 actual=${actual}"
        exit 1
    fi
}

test_all_revc() {
    test_revc '' ''
    test_revc 'A' 'T'
    test_revc 'C' 'G'
    test_revc 'G' 'C'
    test_revc 'T' 'A'
    test_revc 'AAAACCCGGT' 'ACCGGGTTTT'
    echo "REVC: PASS"
}

# test_all_revc



# fib0() {
#     local n=$1
#     local k=$2
#     ((k < 1)) && k=1
# 
#     local a=1
#     local b=1
#     local temp
#     for ((i = 2; i < n; i++)); do
#         temp=$b
#         b=$((k * a + b))
#         a=$temp
#     done
#     echo "$b"
# }

fib() {
    local n=$1
    local k=$2
    ((k < 1)) && k=1

    local arr=(1 1)
    for ((i = 2; i < n; i++)); do
        arr[i]=$((k * arr[i - 2] + arr[i - 1]))
    done
    echo "${arr[${#arr[@]} - 1]}"
}

test_fib() {
    local actual
    actual=$(fib $1 $2)
    if [ "$actual" != "$3" ]; then
        echo "failed fib=$1 $2 expected=$3 actual=${actual}"
        exit 1
    fi
}

test_all_fib() {
    test_fib 1 1 1 
    test_fib 2 1 1 
    test_fib 3 1 2 
    test_fib 4 1 3 
    test_fib 5 1 5 
    test_fib 1 99 1 
    test_fib 2 99 1 
    test_fib 3 2 3 
    test_fib 3 3 4 
    test_fib 5 3 19
    echo "FIB : PASS"
}

# test_all_fib



# gc3() {
#     local tags=( $(echo "$1" | grep '>' | sed 's/>//') )
#     local pattern='>R[a-z]\{7\}_[0-9]\{4\}'
#     local dnas=( $(echo "${1//[[:space:]]/}" | sed "s/${pattern}/ /g") )
#     local maxi=0
#     local max1=0
#     local max2=0
#     for ((i = 0; i < ${#dnas[@]}; i++)); do
#         max2=$(dna "${dnas[i]}" | awk '{printf "%0.6f", ($2+$3) / ($1+$2+$3+$4) * 100}')
#         [[ "$max1" != "$max2" ]] && max1="$max2" && maxi="$i"
#     done
#     echo "${tags[$maxi]}"
#     echo "$max2"
# }

gc() {
    local tags=( $(grep '>' <<< "$1") ) # Alternatives?

    mapfile -t arr <<< "$1"
    local dnas
    dnas="${arr[*]/>*/>}"
    dnas="${dnas// /}"
    IFS=" " read -ra dnas <<< "${dnas//>/ }"

    local v=0
    local maxi=0
    local maxv=0
    local greater
    for ((i = 0; i < ${#dnas[@]}; i++)); do
        v=$(dna "${dnas[i]}" | awk '{printf "%0.6f", ($2+$3) / ($1+$2+$3+$4) * 100}')
        greater=$(echo "$maxv $v" | awk '{printf "%s", $1 < $2 ? 1 : 0}')
        ((greater)) && maxv=$v && maxi=$i
    done
    echo "${tags[$maxi]#*>}"
    echo "$maxv"
}

# gc ">Rosalind_6404
# CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
# TCCCACTAATAATTCTGAGG
# >Rosalind_5959
# CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
# ATATCCATTTGTCAGCAGACACGC
# >Rosalind_0808
# CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
# TGGGAACCTGCGGGCAGTAGGTGGAAT"
# printf "Rosalind_0808\n60.919540\n"



hamm() {
    mapfile -t ss <<< "$1"
    local dist=0
    for ((i = 0; i < ${#ss[0]}; i++)); do
        [[ "${ss[0]:$i:1}" != "${ss[1]:$i:1}" ]] && ((++dist))
    done
    echo "$dist"
}

# hamm "GAGCCTACTAACGGGAT
# CATCGTAATGACGGCCT"
# echo "7"



iprb() {
    local total=$(( "${1// /+}" ))
    local num=$(echo "$1 $total" | awk '{print 4*$1*($2+$3+$4-1) + $2*(4*$3+3*$2-3)}')
    local denom=$((4 * total * (total - 1)))
    echo "$num $denom" | awk '{printf "%0.5f\n", $1 / $2}'
}

# iprb "2 2 2"
# echo "0.78333"



prot() {
    declare -A table=(
        ["UUU"]="F" ["CUU"]="L" ["AUU"]="I" ["GUU"]="V"
        ["UUC"]="F" ["CUC"]="L" ["AUC"]="I" ["GUC"]="V"
        ["UUA"]="L" ["CUA"]="L" ["AUA"]="I" ["GUA"]="V"
        ["UUG"]="L" ["CUG"]="L" ["AUG"]="M" ["GUG"]="V"
        ["UCU"]="S" ["CCU"]="P" ["ACU"]="T" ["GCU"]="A"
        ["UCC"]="S" ["CCC"]="P" ["ACC"]="T" ["GCC"]="A"
        ["UCA"]="S" ["CCA"]="P" ["ACA"]="T" ["GCA"]="A"
        ["UCG"]="S" ["CCG"]="P" ["ACG"]="T" ["GCG"]="A"
        ["UAU"]="Y" ["CAU"]="H" ["AAU"]="N" ["GAU"]="D"
        ["UAC"]="Y" ["CAC"]="H" ["AAC"]="N" ["GAC"]="D"
        ["UAA"]="" ["CAA"]="Q" ["AAA"]="K" ["GAA"]="E"
        ["UAG"]="" ["CAG"]="Q" ["AAG"]="K" ["GAG"]="E"
        ["UGU"]="C" ["CGU"]="R" ["AGU"]="S" ["GGU"]="G"
        ["UGC"]="C" ["CGC"]="R" ["AGC"]="S" ["GGC"]="G"
        ["UGA"]="" ["CGA"]="R" ["AGA"]="R" ["GGA"]="G"
        ["UGG"]="W" ["CGG"]="R" ["AGG"]="R" ["GGG"]="G"
    )

    echo "$1" | fold -w3 | while read -r codon; do printf '%s' "${table[$codon]}"; done
    printf '\n'
}

# prot "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
# echo "MAMAPRTEINSTRING"



subs() {
    local arr=( $(echo "$1") )
    local str=( $(echo "${arr[0]}" | sed 's/./& /g') )
    local substr=( $(echo "${arr[1]}" | sed 's/./& /g') )
    for ((i = 0; i < $(( ${#str[@]} - ${#substr[@]} + 1 )); i++)); do
        local match=0
        for ((j = 0; j < ${#substr[@]}; j++)); do
            [[ "${substr[j]}" != "${str[i + j]}" ]] && break
            ((++match))
        done
        ((match == 4)) && printf "$((++i)) "
    done
    printf '\n'
}

# subs "GATATATGCATATACTT
# ATAT"
# echo "2 4 10"



cons() {
    local seqs=( $(echo "$1" | grep -v '>') )
    local n=${#seqs[@]}
    local m=${#seqs[0]}

    local aa=( $(for ((i=0; i<m; i++)); do echo 0; done) )
    local cc=( $(for ((i=0; i<m; i++)); do echo 0; done) )
    local gg=( $(for ((i=0; i<m; i++)); do echo 0; done) )
    local tt=( $(for ((i=0; i<m; i++)); do echo 0; done) )

    for ((i = 0; i < n; i++)); do
        local seq=( $(echo "${seqs[i]}" | sed 's/./& /g') )
        for ((j = 0; j < m; j++)); do
            [[ "${seq[j]}" == 'A' ]] && ((++aa[j])) && continue
            [[ "${seq[j]}" == 'C' ]] && ((++cc[j])) && continue
            [[ "${seq[j]}" == 'G' ]] && ((++gg[j])) && continue
            [[ "${seq[j]}" == 'T' ]] && ((++tt[j])) && continue
        done
    done

    for ((i = 0; i < ${#aa[@]}; i++)); do
        ((aa[i] >= cc[i] && aa[i] >= gg[i] && aa[i] >= tt[i])) && printf 'A' && continue
        ((cc[i] >= gg[i] && cc[i] >= tt[i])) && printf 'C' && continue
        ((gg[i] >= tt[i])) && printf 'G' && continue
        printf 'T'
    done
    printf '\n'
    echo "A: ${aa[@]}"
    echo "C: ${cc[@]}"
    echo "G: ${gg[@]}"
    echo "T: ${tt[@]}"
}

# cons ">Rosalind_1
# ATCCAGCT
# >Rosalind_2
# GGGCAACT
# >Rosalind_3
# ATGGATCT
# >Rosalind_4
# AAGCAACC
# >Rosalind_5
# TTGGAACT
# >Rosalind_6
# ATGCCATT
# >Rosalind_7
# ATGGCACT"
# echo "ATGCAACT
# A: 5 1 0 0 5 5 0 0
# C: 0 0 1 4 2 0 6 1
# G: 1 1 6 3 0 1 0 0
# T: 1 5 0 0 0 1 1 6"



fibd() {
    local n=$1
    ((n < 2)) && echo 1 && return 0
    local m=$2
    ((m < 1)) && m=1
    ((m < 3)) && echo $((m - 1)) && return 0

    local arr=(1 1)
    for ((i = 2; i < m; i++)); do
        arr[i]=$((arr[i - 1] + arr[i - 2]))
    done
    for ((i = m; i < n; i++)); do
        arr[i]=0
        for ((j = 0; j < m - 1; j++)); do
            arr[i]=$((arr[i] + arr[i - 2 - j]))
        done
    done
    # echo "${arr[${#arr[@]} - 1]}"
    echo "${arr[$n - 1]}"
}

test_fibd() {
    local actual
    actual=$(fibd "$1" "$2")
    if ((actual != $3)); then
        echo "failed n=$1 m=$2 expected=$3 actual=${actual}"
        exit 1
    fi
}

test_all_fibd() {
    test_fibd 1 1 1

    for i in {2..10}; do
        test_fibd 1 "$i" 1
        test_fibd "$i" 1 0
        test_fibd "$i" 2 1
    done

    test_fibd 2 3 1
    test_fibd 3 3 2
    test_fibd 4 3 2
    test_fibd 5 3 3
    test_fibd 6 3 4

    test_fibd 6 4 6

    test_fibd 6 5 7

    echo "FIBD: PASS"
}

# test_all_fibd



grph() {
    local k=3

    local tags=( $(grep '>' <<< "$1" | sed 's/>//g') )

    mapfile -t arr <<< "$1"
    local dnas
    dnas="${arr[*]/>*/>}"
    IFS=" " read -ra dnas <<< "${dnas//>/ }"

    local pref
    local suff
    local n=0
    for ((i = 0; i < ${#dnas[@]}; i++)); do
        suff="${dnas[i]:${#dnas[i]}-$k}"
        for ((j = i + 1; j < ${#dnas[@]}; j++)); do
            pref="${dnas[j]:0:$k}"
            [[ "$suff" == "$pref" ]] && echo "${tags[i]} ${tags[j]}" && ((++n))
        done
    done
}

# grph ">Rosalind_0498
# AAATAAA
# >Rosalind_2391
# AAATTTT
# >Rosalind_2323
# TTTTCCC
# >Rosalind_0442
# AAATCCC
# >Rosalind_5013
# GGGTGGG"
# echo "Rosalind_0498 Rosalind_2391
# Rosalind_0498 Rosalind_0442
# Rosalind_2391 Rosalind_2323"



iev() {
    local AA_AA="$1" # 4/4
    local AA_Aa="$2" # 4/4
    local AA_aa="$3" # 4/4
    local Aa_Aa="$4" # 3/4
    local Aa_aa="$5" # 2/4
    local aa_aa="$6" # 0/4

    local n=2
    local exp_x4=$(((AA_AA + AA_Aa + AA_aa) * 4 + Aa_Aa * 3 + Aa_aa * 2))
    echo $((exp_x4 * n)) | awk '{print $1 / 4}'
}

# iev 1 0 0 1 0 1
# echo 3.5



lcs() {
    declare -A arr
    local z=0
    local res

    for ((i = 0; i < ${#1}; i++)); do
        for ((j = i; j < ${#2}; j++)); do
            if [[ "${1:i:1}" == "${2:j:1}" ]]; then
                if ((i == 0 || j == 0)); then
                    arr[$i,$j]=1
                else
                    local ii=$((i - 1))
                    local jj=$((j - 1))
                    arr[$i,$j]=$((arr[$ii,$jj] + 1))
                fi

                if ((z < ${arr[$i,$j]})); then
                    z=${arr[$i,$j]}
                    res="${1:$((i-z+1)):1}"
                    # echo "$res i=$i j=$j"
                elif ((z == ${arr[$i,$j]})); then
                    res="$res""${1:$((i-z+1)):1}"
                fi
                echo "i=$i j=$j ${1:i:1} z=$z res=$res"
            else
                arr[$i,$j]=0
            fi
            # echo "${arr[@]}"
            # echo "${!arr[@]}"
        done
    done
    echo $z
    echo "$res"
}

lcsm() {
    local x="${1//$'\n'/ }"
    declare -A arr
    declare $(awk 'BEGIN{RS=">"} {print "arr["$1"]="$2}' <<< "${x#>}")

    for x in "${!arr[@]}"; do
        echo "$x: ${arr[$x]}"
    done

    lcs "${arr[Rosalind_1]}" "${arr[Rosalind_2]}"
}

# lcsm ">Rosalind_1
# GATTACA
# >Rosalind_2
# TAGACCA
# >Rosalind_3
# ATACA"
# echo AC



exit



