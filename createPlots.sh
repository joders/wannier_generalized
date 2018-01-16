#!/bin/bash

#gnuplot -e "set terminal 'wxt' font 'Helvetica,16'; unset key; set xlabel 'bandindex'; set ylabel 'overlap'; plot 'overlapWithUlf' pt 7 pointsize 2;" -

#userPlotOptions="with lines"
#userPlotOptions="pointsize 2"
##usetPlotOptions="pt 7 pointsize 2"

#userCommands="set logscale cb; set cbrange [0.01:10];"

#                set xlabel "eigenstate index" font ",20" offset 0,-1 ;
#                set xtics font ",15" ;
#                set ylabel "band index" font ",20" offset -2,0 ;
#                set ytics font ",15" ;
#                set cblabel "band projector expectation" font ",20" offset 4,0 ;
#                set cbtics font ",15" ;
#                set size .9,.86 ;
#                set origin 0,.03 ;
#                set logscale cb ;

for c in {0..89}; do
#for c in {0..59}; do
#for c in {0..89}; do
#for c in 0 1 5 7; do
    dataSetsToPlot[$c]=$c
done

postscriptOptions="landscape enhanced color dashed \"Helvetica\" 24"
#postscriptOptions="landscape enhanced color solid \"Helvetica\" 24"

gpfMargin=( 0.15 .98 0.05 1 )
gpmMargin=( .1 .78 0.05 1 )
gpMargin=( 0.05 .98 0 1 )
margin=( 0 1 0 1 )

makeSetMarginCommands() {
    setMarginCommands="
        set lmargin at screen $1;
        set rmargin at screen $2;
        set bmargin at screen $3;
        set tmargin at screen $4;
        "
}

makeGnuplotCommandList() {  # 1.arg: rel. path to datafile
    if test "`sed -n '1p' "$1" | cut -f1`" == "#"; then 
        title=`sed -n '1p' "$1" | cut -f2`
        xaxis=`sed -n '1p' "$1" | cut -f3`
        yaxis=`sed -n '1p' "$1" | cut -f4`
        cbaxis=`sed -n '1p' "$1" | cut -f5`
        if test -z "$userPlotOptions"; then
            plotOptions="`sed -n '1p' "$1" | cut -f6`"
        else
            plotOptions="$userPlotOptions"
        fi
        commandsFromFile=`sed -n '1p' "$1" | cut -f7`
#        echo $title
#        echo $xaxis
#        echo $yaxis
        setLabelsCommands="
            set title \"$title\";
            set xlabel \"$xaxis\";
            set ylabel \"$yaxis\";
            set cblabel \"$cbaxis\";";
    else
        setLabelsCommands="";
    fi

    filepath="$1"
    filename=$(basename "$1")
    extension="${filename##*.}"
    filename="${filename%.*}"
    if test "$extension" = "$filename"; then
        extension=""
    fi
#    echo $filename
#    echo $extension

    setKeyCommand="unset key;"

    if test "$extension" = "gp"; then
        plotArguments="\"$1\" $plotOptions,"
        makeSetMarginCommands ${gpMargin[@]}
    else
        if test "$extension" = "gpf"; then
            makeSetMarginCommands ${gpfMargin[@]}
            if test ${#scriptArguments[@]} -ge 2; then
                dataSetsToPlot=("${scriptArguments[@]:1}")
            fi

            if test `sed -n '2p' "$1" | cut -f1` == "#"; then 
                setKeyCommand="set key box opaque;"
                (( index=2 ))
                while true; do
                    nextLabel=`sed -n '2p' "$1" | cut -f$index`
                    if test -z "$nextLabel"; then
                        break
                    fi
                    setLegendLabels[ (($index-2)) ]=$nextLabel
                    (( index=$index+1 ))
                done
                plotArguments=""
                for c in "${dataSetsToPlot[@]}"; do
                    plotArguments="$plotArguments
                                   \"$1\" using 1:$(( $c+2 )) title \"${setLegendLabels[$c]}\" $plotOptions,"
                done
            else
                plotArguments=""
                for c in "${dataSetsToPlot[@]}"; do
                    plotArguments="$plotArguments
                                   \"$1\" using 1:$(( $c+2 )) $plotOptions,"
                done
            fi
        else 
            if test "$extension" = "gpm"; then
                makeSetMarginCommands ${gpmMargin[@]}

                if [[ $commandsFromFile =~ .*logscale.* ]] || [[ $userCommands =~ .*logscale.* ]]; then
                    filepath="${filepath}.tmpabs"
                    vim -c':2,$ s/ -/ /g' -c"w $filepath" -c'q!' "$1"
                fi

                if test ${#scriptArguments[@]} -eq 3; then
                    fI=${scriptArguments[1]} #firstIndex
                    lI=${scriptArguments[2]} #lastIndex
                    plotArguments="\"$filepath\" matrix every ::$fI:$fI:$lI:$lI with image $plotOptions"
                else
                    plotArguments="\"$filepath\" matrix with image $plotOptions"
                fi
            fi
        fi
    fi

    setTerminalCommand="
        set term postscript $postscriptOptions;
        set output \"$filepath.ps\";"

    gnuplotCommandList="
        $commandsFromFile
        $setLabelsCommands
        $setKeyCommand
        $setMarginCommands
        $userCommands
        $setTerminalCommand
        plot $plotArguments"

    #echo $gnuplotCommandList
}


# START #

scriptArguments=("$@")

if test "$1" = ""; then
    plotDir="output"
else
    if test -d "$1"; then
        plotDir="$1"
    else 
        plotDir=""
        if test -f "$1"; then
            makeGnuplotCommandList "$1"
            gnuplot -e "$gnuplotCommandList" 2> /dev/zero
            mv "${filepath}.ps" "${filepath}.ps_wrongColor"
            patch -s -o "${filepath}.ps" "${filepath}.ps_wrongColor" ../psColor.patch
            trash "${filepath}.ps_wrongColor"
            ps2pdf "${filepath}.ps"
            trash "${filepath}.ps"
            /usr/bin/okular `basename "${filepath}.pdf"`>/dev/zero 2>/dev/zero
            trash `basename "${filepath}.pdf"`
        else
            echo "something strange happened"
        fi
    fi
fi

if test "$plotDir" != ""; then
    cd "$plotDir"
    for i in *.gp*; do
        if test -e "$i"; then
            makeGnuplotCommandList "$i"
            gnuplot -e "$gnuplotCommandList" 2> /dev/zero
            mv "$i.ps" "$i.ps_wrongColor"
            patch -s -o "$i.ps" "$i.ps_wrongColor" ../psColor.patch
            trash "$i.ps_wrongColor"
            ps2pdf "$i.ps"
            trash "$i.ps"
        fi
    done
fi
