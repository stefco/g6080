rm -f real_time_clock_o3_stef_assignment_1
rm -f real_time_clock_o0_stef_assignment_1
rm -f timing-results-o3.csv
rm -f timing-results-o0.csv

gcc -O3 -o q2_real_time_clock_o3_stef_assignment_1 real_time_clock_stef_assignment_1.c
gcc -O0 -o q2_real_time_clock_o0_stef_assignment_1 real_time_clock_stef_assignment_1.c

echo doing twenty loops

for i in `seq 1 20`; do
	./q2_real_time_clock_o3_stef_assignment_1 >> timing-results-o3.csv
	./q2_real_time_clock_o0_stef_assignment_1 >> timing-results-o0.csv
done
