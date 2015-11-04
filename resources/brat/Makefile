all: brat brat-large trim acgt-count brat-large-build rev-compl check-strands remove-dupl convert-to-sam
brat: BRAT.cpp
	g++ -O3 -w -o brat BRAT.cpp
brat-large: BRAT_large.cpp
	g++ -O3 -w -o brat-large BRAT_large.cpp
brat-large-build: BRAT_large_build.cpp
	g++ -O3 -w -o brat-large-build BRAT_large_build.cpp
trim: Filter_and_trim.cpp
	g++ -O3 -w -o trim Filter_and_trim.cpp
acgt-count: Coverage_ACGT.cpp
	g++ -O3 -w -o acgt-count Coverage_ACGT.cpp
rev-compl: Take_reverse_complement.cpp
	g++ -O3 -w -o rev-compl Take_reverse_complement.cpp
check-strands: Check_strands.cpp
	g++ -O3 -w -o check-strands Check_strands.cpp
remove-dupl: Remove_copy_duplicates.cpp
	g++ -O3 -w -o remove-dupl Remove_copy_duplicates.cpp
convert-to-sam: Convert_format_to_SAM.cpp
	g++ -O3 -w -o convert-to-sam Convert_format_to_SAM.cpp

clean: brat brat-large trim acgt-count brat-large-build rev-compl check-strands remove-dupl convert-to-sam
	rm brat brat-large trim acgt-count brat-large-build rev-compl check-strands remove-dupl convert-to-sam
