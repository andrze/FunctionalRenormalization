
countoutput=`ls -1 *.sh.o* 2>/dev/null | wc -l`
counterror=`ls -1 *.sh.e* 2>/dev/null | wc -l`
if [ $(($countoutput + $counterror)) -gt 0 ]
then
    echo "Cleaning log files"
    rm *.sh.o* *.sh.e*
else
    echo "No log files to clean"
fi


if [ $# -gt 0 ]
then 
    if [ "${1: -1}" == "/" ]
    then
        directory="${1:0: ${#1}-1}"
        
    else
        directory="$1"
    fi

    if [ -d $directory ]
    then
        
        echo "Cleaning the '$directory' directory"
        rm -r "${directory}"
        
    else
        echo "The '$directory' directory does not exist"
    fi
fi 


    
