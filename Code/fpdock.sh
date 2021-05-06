set -e
#<<'COMMENTS'
#recc=['A']
#ligc=['B']
#ja="4m76_AB"
lis=(Data/[1-9]*)
for ja in "${lis[@]}"
do
	cry="$(cut -d'/' -f2 <<<"$ja").pdb"  #"${ja}.pdb" #"${ja}.pdb"
	sub="$(cut -d'_' -f1 <<<$(cut -d'/' -f2 <<<"$ja"))"  #${ja:$((${#ja}-7)):4}
	mkdir ${sub}
for h in {0..0}
do
	#cp "$j/${j:$((${#j}-7)):4}_u1.pdb" ./
	#cp "$j/${j:$((${#j}-7)):4}_u2.pdb" ./
	#cp "$j/${j:$((${#j}-7)):7}_st.pdb" ./
	#cp "$j/${j}_u1.pdb" ./
	#cp "$j/${j}_u2.pdb" ./
	#cp "$j/${j}.pdb" ./

	
	touch  "report_${sub}_${h}.txt"
	#bash -c '[ -d ./poses ] && rm -r ./poses/*'
	#python 30-6-20-fpa.py -rpdb "${ja}/${ja}_u1.pdb" -lpdb "${ja}/${ja}_u2.pdb" -recchain ${recc} -ligchain ${ligc} -n 30
	chain=($(python3 dock.py -pdb ${ja} -n 200 | tr -d '[],'))
	echo ${chain[0]}

	ch=()
	for ci in $(seq 1 ${chain[0]})
	do
		eval a=${chain[ci]}
		ch+=(${a})
	done
	echo ${ch[@]}

	result=(native_${sub}/*)

	# ${ch[1]} ${ch[2]} ${ch[3]} ${ch[4]} ${ch[5]}

	for i in "${result[@]}"
	do
		./DockQ.py ${ja}/${cry} $i  -native_chain1 ${ch[@]}  | tee -a ds.txt # -native_chain1 -native_chain1 A B
		result1=($(python read.py | tr -d '[],'))
		result1=("${result1[@]//\'/}")
		echo "$i  ${result1[@]}" >> "report_${sub}_${h}.txt"
		#printf "\n" >>"result_${sub}.txt"
		#arr += ${result[@]}

		rm ds.txt
	done

	tar cvfz "native_${sub}_${h}.tar.gz" "native_${sub}"
	mv native_${sub}_${h}.tar.gz report_${sub}_${h}.txt "${sub}"
	
	rm -r poses/*
	rm -r native_${sub}
	rm -r /tmp/tmp[a-d]*
	rm -r /tmp/tmp[e-j]*
	rm -r /tmp/tmp[k-r]*
	rm -r /tmp/tmp[s-z]*
	rm -r /tmp/tmp[0-9]*
	rm -r /tmp/tmp*
	#rm -r tmp*

done

touch "best_${sub}_${h}.txt"

li=(${sub}/report*)
for xi in "${li[@]}"
do
	mapfile -t array < <(python final.py -f ${xi})
	echo ${array[0]} >> "best_${sub}.txt"
done
mv best_${sub}.txt bestenergy.txt "${sub}"
rm -r *.pdb
rm -r Data/${sub}*
#rm -r /tmp/tmp*
#mv ${sub} "/home/sharon/IFPAresult50"
#mv ${ja} "/home/sharon/done_PPDBv5"
done


