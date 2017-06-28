#!/usr/bin/env bash
set -u

exit_status=0

smof_src=https://github.com/arendsee/smof
synder_src=https://github.com/arendsee/synder

pconf=preconfig.sh
rconf=runconfig.R

[[ -d bin ]] || mkdir bin

make-config(){
    base=$1
    dest=$2
    echo -n "Checking for $base ... "
    tconf=etc/template-$base
    if [[ ! -r $base ]]
    then
        echo "initializing"
        sed -n '/^$/,$ p' $tconf | sed "s;FAGIN_HOME;$PWD;" > $base
        ln -sf $PWD/$base $dest
    else
        echo "already exists"
    fi
}

check-exe(){
    echo -n "Looking for $1 ... "
    if [[ -z `type -P $1` ]]
    then
        echo "MISSING"
        echo $2
        exit_status=1
    else
        echo "OK"
    fi
}

install-exe(){
    echo -n "Looking for $1 ... "
    inpath=0 inbin=0

    if [[ ! -z `type -P $1` ]]
    then
        inpath=1
        path=`which $1`
    fi

    if [[ -f "bin/$1" ]]
    then
        inbin=1
        path="$PWD/bin/$1"
    fi

    if [[ $inpath -eq 1 || $inbin -eq 1 ]]
    then
        echo "OK"
        sed "s;^$1=.*;$1='$path';" $pconf > $pconf~
        mv $pconf~ $pconf
    else
        echo "MISSING"
        $2 $1
        if [[ $? -ne 0  ]]
        then
            echo "Failed to install $1"
            exit_status=1
        fi
    fi
}

install-smof(){
    if [[ -f bin/smof ]]
    then
        return 0
    else
        d=src/smof
        git clone $smof_src $d || return 1
        cp $d/smof.py bin/smof
        rm -rf $d
    fi
    return 0
}

install-synder(){
    d=src/synder
    if [[ -f bin/synder ]]
    then
        return 0
    else
        git clone -b v0.17.4 $synder_src $d || return 1
        cd $d
        make || return 1
        cd -
        cp $d/synder bin
        [[ -d $d ]] && rm -rf $d
    fi
    return 0
}

make-config $pconf $PWD/src/prologue/config
make-config $rconf $PWD/src/report/config

check-exe python3  "python3 not in path, please install"
check-exe latexmk  "latexmk not in path (needed for report generation), please install"
check-exe R        "R not in path, please install"
check-exe parallel "parallel is not in your path"
check-exe Rscript  "Rscript not in path, please install R"
check-exe bedtools "bedtools not in path, please install"
check-exe transeq  "emboss::transeq not in path, please install emboss"
check-exe getorf   "emboss::getorf not in path, please install emboss"
check-exe git      "git not in path, please install"

install-exe smof   install-smof
install-exe synder install-synder

./src/prologue/configure.R

if [[ $exit_status -eq 0 ]]
then
cat << EOF
 Configuration successful, you may now:
   1. Set required fields in preconfig.sh
   2. Set required fields in runconfig.R
   3. make load  # takes several minures
   4. make run   # this may take a few hours
EOF
else
    echo "!!! Configuration FAIL: do NOT proceed, do NOT run make !!!"
fi

exit $exit_status
