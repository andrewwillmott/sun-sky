

mkdir -p images

SUNSKY="./sunsky -f -v -l 30 0 -d 58 -b 3"
OPTS=

TIMES="12.5 13.5 14.5 15.5 16.5 17.5 18.5"
ROUGHNESSES="0 0.166 0.333 0.5 0.666 0.833 1"

for type in hosek preetham; do
    echo $type ":"
    i=1
    for time in $TIMES; do
        $SUNSKY $OPTS -s $type -t $time && mv -f sky-hemi.png images/${type}-${i}.png
        let "i++"
    done
done

OPTS="-o 0.5"

for type in hosek preetham; do
    echo $type ":"
    i=1
    for time in $TIMES; do
        $SUNSKY $OPTS -s $type -t $time && mv -f sky-hemi.png images/${type}-oc-${i}.png
        let "i++"
    done
done

OPTS=""

for type in hosek preetham; do
    echo $type ":"
    i=1
    for o in $ROUGHNESSES; do
        $SUNSKY $OPTS -s $type -t 14 -o $o && mv -f sky-hemi.png images/${type}-ocd-${i}.png
        # $SUNSKY $OPTS -s $type -t 07 -o $o && mv -f sky-hemi.png images/${type}-ocs-${i}.png
        let "i++"
    done
done

OPTS=""

for type in hosekBRDF preethamBRDF; do
    echo $type ":"
    i=1
    for r in $ROUGHNESSES; do
        $SUNSKY $OPTS -s $type -t 14 -r $r && mv -f sky-hemi.png images/${type}-rd-${i}.png
        $SUNSKY $OPTS -s $type -t 07 -r $r && mv -f sky-hemi.png images/${type}-rs-${i}.png
        let "i++"
    done
done

mogrify -resize 120x120 images/*
