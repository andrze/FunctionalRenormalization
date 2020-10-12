rm *.sh.o*
rm *.sh.e*

if [ $# -gt 0 ]
then 
    if [ "${1: -1}" == "/" ]
    then
        directory="${1::-1}"
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

