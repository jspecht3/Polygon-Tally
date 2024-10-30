if [ $# -eq 0 ];
	then
		echo "You need 4 parameters: sides of polygon, polygon radius, x coordinate, y coordinate."
		exit 1

elif [ $# -gt 4 ];
	then
		echo "Only 4 paramters needed: sides of polygon, polygon radius, x coordinate, y coordinate."
		exit 1

else
	n=$1
	rad=$2
	x=$3
	y=$4

	echo ""
	echo "----- Inputs -----"
	echo "Number of Sides...: $n"
	echo "Polygon Radius....: $rad"
	echo "x coordinate......: $x"
	echo "y coordinate......: $y"
fi

echo "$n,$rad,$x,$y" >> temp-check.txt

python r-alpha.py


mv cpp-implementation/main .
echo "$n $rad $x $y" | ./main

mv main cpp-implementation/

rm temp-check.txt
