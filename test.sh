for X in 10 100 1000 10000 100000 1000000 10000000
do  
    for R in 1 2 4 8 16
    do
        echo X = $X, R = $R
        ./db5242 $X 0 0 0 $R
        echo
    done
done