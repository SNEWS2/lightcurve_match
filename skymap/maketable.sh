getval(){
    echo `tail -n 17 $1 | head -n $2 | tail -n 1 | cut -f $3 -d ' '`
}
fillrow() {
    val1=`getval log_256_100000_GC_IC_HK_JUNO.txt $1 $2`
    val1=`printf %.1f $val1`
    val2=`getval log_256_100000_GC_IC_HK_ARCA.txt $1 $2`
    val2=`printf %.1f $val2`
    val3=`getval log_256_100000_GC_IC_ARCA_JUNO.txt $1 $2`
    val3=`printf %.1f $val3`
    val4=`getval log_256_100000_GC_HK_ARCA_JUNO.txt $1 $2`
    val4=`printf %.1f $val4`
    val5=`getval log_256_100000_GC_IC_HK_ARCA_JUNO.txt $1 $2`
    val5=`printf %.1f $val5`
    echo "$val1 & $val2 & $val3 & $val4 & $val5"
}
fillrowerr() {
    val1=`getval log_256_100000_GC_IC_HK_JUNO.txt $1 $2`
    val1=`printf %.1f $val1`
    err1=`getval log_256_100000_GC_IC_HK_JUNO.txt $1 $3`
    err1=`printf %.1f $err1`

    val2=`getval log_256_100000_GC_IC_HK_ARCA.txt $1 $2`
    val2=`printf %.1f $val2`
    err2=`getval log_256_100000_GC_IC_HK_ARCA.txt $1 $3`
    err2=`printf %.1f $err2`

    val3=`getval log_256_100000_GC_IC_ARCA_JUNO.txt $1 $2`
    val3=`printf %.1f $val3`
    err3=`getval log_256_100000_GC_IC_ARCA_JUNO.txt $1 $3`
    err3=`printf %.1f $err3`

    val4=`getval log_256_100000_GC_HK_ARCA_JUNO.txt $1 $2`
    val4=`printf %.1f $val4`
    err4=`getval log_256_100000_GC_HK_ARCA_JUNO.txt $1 $3`
    err4=`printf %.1f $err4`

    val5=`getval log_256_100000_GC_IC_HK_ARCA_JUNO.txt $1 $2`
    val5=`printf %.1f $val5`
    err5=`getval log_256_100000_GC_IC_HK_ARCA_JUNO.txt $1 $3`
    err5=`printf %.1f $err5`

    echo $val1'$\pm$'$err1 \& $val2'$\pm$'$err2 \& $val3'$\pm$'$err3 \& $val4'$\pm$'$err4 \& $val5'$\pm$'$err5
}
echo 'GC ndof=2 \#1 90\% CL area (deg$^2$) &' `fillrow 6 9` '\\\\ \hline'
echo 'GC ndof=2 \#2 90\% CL area (deg$^2$) &' `fillrowerr 2 6 8` '\\\\ \hline'
echo 'GC \#3 90\% C.L. area (deg$^2$) &' `fillrow 4 9` '\\\\ \hline' 
echo 'GC ndof=2 \#1 68\% CL area (deg$^2$) &' `fillrow 5 9` '\\\\ \hline'
echo 'GC ndof=2 \#2 68\% CL area (deg$^2$) &' `fillrowerr 1 6 8` '\\\\ \hline'
echo 'GC \#3 68\% C.L. area (deg$^2$) &' `fillrow 3 9` '\\\\ \hline'

echo 'GC \#2 90\% real coverage (\%) &' `fillrowerr 2 12 14` '\\\\ \hline'
echo 'GC \#2 68\% real coverage (\%) &' `fillrowerr 1 12 14` '\\\\ \hline'
