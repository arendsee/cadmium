print-warning(){
    if [[ -t 2 ]]
    then
        echo -e "\e[1;31mERROR: \e[0m$1"
    else
        echo $1
    fi
}

croak(){
    print-warning "$1"
    exit 1
}

check-read(){
    file="$1"   # the name of the test file
    place="$2"  # file context, e.g. foo.sh::doit (for debugging)
    if [[ ! -r "$file" ]]
    then
        croak "Cannot read '$file', used in $place"
    fi
}

check-dir(){
    dir=$1    # the name of the test file
    place=$2  # dir context, e.g. foo.sh::doit (for debugging)
    if [[ ! -r "$dir" ]]
    then
        croak "Cannot access directory '$dir', used in $place"
    fi
}

safe-mkdir(){
    dir=$1
    [[ -d $dir ]] || mkdir -p $dir || croak "Could not make directory '$dir'"
}

kill-mkdir(){
    dir=$1
    [[ -d $dir ]] && rm -rf $dir
    mkdir $dir
}

# `bedtools getfasta` reads the sequence name from the third column of the gff.
# This function moves some other column (ARG1) to the 3rd.
rename_for_bedtools (){
    awk -v from=$1 '
        BEGIN{FS="\t"; OFS="\t"}
        {$3 = $from}
        {print}
    '
}
