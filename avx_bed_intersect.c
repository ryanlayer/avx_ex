#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <immintrin.h>

#define GENOME_SIZE 3095677412

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

static struct timeval _start, _stop;

void start()
{
    gettimeofday(&_start,0);
}

void stop()
{
    gettimeofday(&_stop,0);
}

unsigned long report()
{
    return (_stop.tv_sec - _start.tv_sec) * 1000000 + //sec to microsec
            _stop.tv_usec - _start.tv_usec;
}

uint32_t bed_to_bits(char *file_name,
                     uint32_t **B,
                     uint32_t G_size)
{ 

    uint32_t ints_per_genome = (G_size + 32 - 1)/32;
    int b = posix_memalign((void **)B,
                            32,
                            ints_per_genome*sizeof(unsigned int));
    memset(*B, 0, ints_per_genome*sizeof(unsigned int));

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *pch;

    FILE *bed_file = fopen(file_name, "r");
    if (!bed_file) {
        fprintf(stderr, "Unable to open %s\n", file_name);
        return 1;
    }
 
    uint32_t start,
             end,
             interval_len,
             start_int_i,
             start_bit_i,
             new_bits,
             fill_bits,
             num_intervals = 0;

    while ((read = getline(&line, &len, bed_file)) != -1) {
        if (line[0] != '#') {
            pch = strtok(line, "\t");

            if (pch[0] != '1') {
                fprintf(stderr, "Sorry, this example only cosiders chrom 1");
                return 1;
            }

            // use 0-based closed interval, as it should be
            pch = strtok(NULL, "\t");
            start = atoi(pch);

            pch = strtok(NULL, "\t");
            end = atoi(pch) - 1;  

            start_int_i = start / 32;
            start_bit_i = start - 32*start_int_i;

            interval_len = end - start + 1;

            while (interval_len > 0) {
                fill_bits = MIN(interval_len, (32 - start_bit_i));
                new_bits = (0xffffffff >> (32 - fill_bits)) 
                            << (32 - start_bit_i - fill_bits);
                (*B)[start_int_i] |= new_bits;

                start_int_i += 1;
                interval_len -= fill_bits;
                start_bit_i = 0;
            }
            num_intervals += 1;
        }
    }
    fclose(bed_file);
    free(line);

    return num_intervals;
}

uint32_t intersect_beds(uint32_t *A,
                        uint32_t *B,
                        uint32_t **R,
                        uint32_t G_size)
{
    uint32_t ints_per_genome = (G_size + 32 - 1)/32;

    *R = (uint32_t *) calloc(ints_per_genome, sizeof(uint32_t));

    uint32_t i;

    for (i = 0; i < ints_per_genome; ++i)
        (*R)[i] = A[i] & B[i];

    return ints_per_genome;
}

uint32_t avx_intersect_beds(uint32_t *A,
                            uint32_t *B,
                            uint32_t **R,
                            uint32_t G_size)
{
    uint32_t ints_per_genome = (G_size + 32 - 1)/32;

    *R = (uint32_t *) calloc(ints_per_genome, sizeof(uint32_t));
    int r = posix_memalign((void **)R,
                            32,
                            ints_per_genome*sizeof(unsigned int));

    __m256i *A_avx = (__m256i *)A;
    __m256i *B_avx = (__m256i *)B;
    __m256i *R_avx = (__m256i *)*R;

    __m256 a, b;

    uint32_t i;

    for (i = 0; i < (ints_per_genome/8); ++i) {
        a = _mm256_load_si256(&A_avx[i]);
        b = _mm256_load_si256(&B_avx[i]);
        R_avx[i] = _mm256_and_ps(a,b);
    }

    return ints_per_genome;
}



int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, 
                "usage:\t%s <bed file 1> <bed fille2>\n",
                argv[0]);

        return 1;
    }


    char *bed_1_file_name = argv[1];
    char *bed_2_file_name = argv[2];

    uint32_t *bed_1, *bed_2, *bed_R, *bed_R2;

    uint32_t bed_1_N = bed_to_bits(bed_1_file_name, &bed_1, GENOME_SIZE);
    uint32_t bed_2_N = bed_to_bits(bed_2_file_name, &bed_2, GENOME_SIZE);

    start();
    uint32_t bed_R_N = intersect_beds(bed_1, bed_2, &bed_R, GENOME_SIZE);
    stop();
    printf("%lu\n", report());

    start();
    uint32_t bed_R2_N = avx_intersect_beds(bed_1, bed_2, &bed_R2, GENOME_SIZE);
    stop();
    printf("%lu\n", report());


    uint32_t i;
    for (i = 0; i < bed_R2_N; ++i)
        if (bed_R[i] != bed_R2[i])
            printf("%u\t%u\t%u\n", i, bed_R[i], bed_R2[i]);



    free(bed_1);
    free(bed_2);
    free(bed_R);
    free(bed_R2);
}
