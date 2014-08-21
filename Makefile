all: avx_bed_intersect

avx_bed_intersect: avx_bed_intersect.c
	gcc -mavx -o avx_bed_intersect avx_bed_intersect.c
