I=1
C=200
for ((i=$I;i<$I+$C;++i))
do
#	./dalitz > log_without_LASS/log_$i.txt
	./dalitz > log_floating_smoothEff/log_$i.txt
done
