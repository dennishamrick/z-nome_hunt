/*
Z-HUNT-2 computer program, Sept.19, 1990
*/

/*
serialized i/o to allow to run against large datasets
Updated i/o permissions to match posix. campt 1/10/2000
*/

/*
Turbo C compiler Ver. 1.5
Compact model, 8086/80286 instruction set, math emulator/8087/80287
Speed option, unsigned char
Run on IBM PC XT/AT or compatibles
*/

/*
Written by Ping-jung Chou, under the instruction of Pui S. Ho, according to the
paper "A computer aided thermodynamic approach for predicting the formation of
Z-DNA in naturally occurring sequences",
The EMBO Journal, Vol.5, No.10, pp2737-2744, 1986
With 0.22 kcal/mol/dinuc for mCG (Zacharias et al, Biochemistry, 1988, 2970)
*/

/*Z-NOME-HUNT, October 11 2024
 */
/*
Edited by Dennis Hamrick based off of https://github.com/Ho-Lab-Colostate/zhunt/tree/master/src
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

double linear_search(double x1, double x2, double tole, double (*func)());
double delta_linking(double dl);
double find_delta_linking(int dinucleotides);
double delta_linking_slope(double dl);
void assign_bzenergy_index(int nucleotides, char seq[]);
void best_anti_syn(int dinucleotides, float esum);
void anti_syn_energy(int din, int dinucleotides, float esum);
double assign_probability(double dl);
void calculate_zscore(double a, int maxdinucleotides, int min, int max, char *filename, int startsite, char *chrom);
FILE *open_file(int mode, char *filename, char *typestr);
unsigned input_sequence(FILE *file, int nucleotides, int showfile);

int terms;
double *bztwist, *logcoef, *exponent;
const double _k_rt = -0.2521201;  /* -1100/4363 */
const double sigma = 16.94800353; /* 10/RT */
double deltatwist;
const double explimit = -600.0;
int ss;
double *bzenergy, *best_bzenergy; /* dinucleotides */
float best_esum;                  /* assigned before call to anti_syn_energy() */
char *best_antisyn, *antisyn;     /* nucleotides */
char *tempstr, *sequence;

/* Delta BZ Energy of Dinucleotide */
double dbzed[4][25] = {
    /* AS-AS */
    /* AA    AT    AG    AC    AN    TA    TT    TG    TC    TN    GA    GT    GG    GC    GN    CA    CT    CG    CC    CN    NA    NT     NG    NC    NN */
    {4.40, 6.20, 3.40, 5.20, 3.00, 9.99, 4.40, 1.40, 3.30, 1.90, 9.99, 5.20, 2.40, 4.20, 1.49, 9.99, 3.40, 0.66, 2.40, 0.93, 9.99, 9.99, 9.99, 9.99, 9.99},
    /* SA-SA */
    {4.40, 2.50, 3.30, 1.40, 0.80, 9.99, 4.40, 5.20, 3.40, 1.90, 9.99, 1.40, 2.40, 0.66, -0.71, 9.99, 3.30, 4.20, 2.40, 2.73, 9.99, 9.99, 9.99, 9.99, 9.99},
    /* AS-SA */
    {6.20, 6.20, 5.20, 5.20, 4.60, 9.99, 6.20, 5.20, 5.20, 4.60, 9.99, 5.20, 4.00, 4.00, 2.53, 9.99, 5.20, 4.00, 4.00, 2.53, 9.99, 9.99, 9.99, 9.99, 9.99},
    /* SA-AS */
    {6.20, 6.20, 5.20, 5.20, 4.60, 9.99, 6.20, 5.20, 5.20, 4.60, 9.99, 5.20, 4.00, 4.00, 2.53, 9.99, 5.20, 4.00, 4.00, 2.53, 9.99, 9.99, 9.99, 9.99, 9.99}};

double expdbzed[4][25]; /* exp(-dbzed/rt) */
int *bzindex;           /* dinucleotides */

double linear_search(double x1, double x2, double tole, double (*func)(double))
{
    double f, fmid, dx, xmid, x;

    f = func(x1);
    fmid = func(x2);
    if (f * fmid >= 0.0)
        return x2;
    x = (f < 0.0) ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
    do
    {
        dx *= 0.5;
        xmid = x + dx;
        fmid = func(xmid);
        if (fmid <= 0.0)
            x = xmid;
    } while (fabs(dx) > tole);
    return x;
}

double delta_linking(double dl)
{
    double sump, sumq, z, expmini;
    int i;

    expmini = 0.0;
    for (i = 0; i < terms; i++)
    {
        z = dl - bztwist[i];
        exponent[i] = z = logcoef[i] + _k_rt * z * z;
        if (z < expmini)
            expmini = z;
    }
    expmini = (expmini < explimit) ? explimit - expmini : 0.0;
    sump = sumq = 0.0;
    for (i = 0; i < terms; i++)
    {
        z = exp(exponent[i] + expmini);
        sumq += z;
        sump += bztwist[i] * z;
    }
    sumq += exp(_k_rt * dl * dl + sigma + expmini);
    return deltatwist - sump / sumq;
}

double find_delta_linking(int dinucleotides)
{
    double sum;
    int i, j;

    for (i = 0; i < dinucleotides; i++)
        bzenergy[i] = 1.0;
    for (i = 0; i < dinucleotides; i++)
    {
        sum = 0.0;
        for (j = 0; j < dinucleotides - i; j++)
        {
            bzenergy[j] *= best_bzenergy[i + j];
            sum += bzenergy[j];
        }
        logcoef[i] = log(sum);
    }
    terms = dinucleotides;
    return linear_search(10.0, 50.0, 0.001, delta_linking);
}

double delta_linking_slope(double dl)
{
    double sump, sump1, sumq, sumq1, x, y, z, expmini;
    int i;

    expmini = 0.0;
    for (i = 0; i < terms; i++)
    {
        z = dl - bztwist[i];
        exponent[i] = z = logcoef[i] + _k_rt * z * z;
        if (z < expmini)
            expmini = z;
    }
    expmini = (expmini < explimit) ? explimit - expmini : 0.0;
    sump = sump1 = sumq = sumq1 = 0.0;
    x = 2.0 * _k_rt;
    for (i = 0; i < terms; i++)
    {
        z = dl - bztwist[i];
        y = exp(exponent[i] + expmini);
        sumq += y;
        sump += bztwist[i] * y;
        y *= z * x;
        sumq1 += y;
        sump1 += bztwist[i] * y;
    }
    y = exp(_k_rt * dl * dl + sigma + expmini);
    sumq += y;
    sumq1 += x * dl * y;
    return (sump1 - sump * sumq1 / sumq) / sumq;
} /* slope at delta linking = dl */

void assign_bzenergy_index(int nucleotides, char seq[])
{
    int i, j, idx;
    char c1, c2;

    i = j = 0;
    do
    {
        c1 = seq[i++];
        c2 = seq[i++];
        switch (c1)
        {
        case 'a':
            switch (c2)
            {
            case 'a':
                idx = 0;
                break;
            case 't':
                idx = 1;
                break;
            case 'g':
                idx = 2;
                break;
            case 'c':
                idx = 3;
                break;
            case 'n':
                idx = 4;
            }
            break;
        case 't':
            switch (c2)
            {
            case 'a':
                idx = 5;
                break;
            case 't':
                idx = 6;
                break;
            case 'g':
                idx = 7;
                break;
            case 'c':
                idx = 8;
                break;
            case 'n':
                idx = 9;
            }
            break;
        case 'g':
            switch (c2)
            {
            case 'a':
                idx = 10;
                break;
            case 't':
                idx = 11;
                break;
            case 'g':
                idx = 12;
                break;
            case 'c':
                idx = 13;
                break;
            case 'n':
                idx = 14;
            }
            break;
        case 'c':
            switch (c2)
            {
            case 'a':
                idx = 15;
                break;
            case 't':
                idx = 16;
                break;
            case 'g':
                idx = 17;
                break;
            case 'c':
                idx = 18;
                break;
            case 'n':
                idx = 19;
            }
            break;
        case 'n':
            switch (c2)
            {
            case 'a':
                idx = 20;
                break;
            case 't':
                idx = 21;
                break;
            case 'g':
                idx = 22;
                break;
            case 'c':
                idx = 23;
                break;
            case 'n':
                idx = 24;
            }
            break;
        }
        bzindex[j++] = idx;
    } while (i < nucleotides);
}

void best_anti_syn(int dinucleotides, float esum)
{
    int i;
    double dl, slope;

    if (esum < best_esum)
    {
        best_esum = esum;
        for (i = 0; i < dinucleotides; i++)
            best_bzenergy[i] = bzenergy[i];
        strcpy(best_antisyn, antisyn);
    }
}

void anti_syn_energy(int din, int dinucleotides, float esum)
{
    int i, nucleotides;
    float e;

    nucleotides = 2 * din;

    antisyn[nucleotides] = 'A';
    antisyn[nucleotides + 1] = 'S';
    i = (din == 0) ? 0 : ((antisyn[nucleotides - 1] == 'S') ? 0 : 3);
    e = dbzed[i][bzindex[din]];
    esum += e;
    bzenergy[din] = expdbzed[i][bzindex[din]];
    if (++din == dinucleotides)
        best_anti_syn(dinucleotides, esum);
    else
        anti_syn_energy(din, dinucleotides, esum);
    esum -= e;
    din--;

    antisyn[nucleotides] = 'S';
    antisyn[nucleotides + 1] = 'A';
    i = (din == 0) ? 1 : ((antisyn[nucleotides - 1] == 'A') ? 1 : 2);
    esum += dbzed[i][bzindex[din]];
    bzenergy[din] = expdbzed[i][bzindex[din]];
    if (++din == dinucleotides)
        best_anti_syn(dinucleotides, esum);
    else
        anti_syn_energy(din, dinucleotides, esum);
}

/* calculate the probability of the value 'dl' in a Gaussian distribution */
/* from "Data Reduction and Error Analysis for the Physical Science" */
/* Philip R. Bevington, 1969, McGraw-Hill, Inc */

double assign_probability(double dl)
{
    static double average = 29.6537135;
    static double stdv = 2.71997;
    static double _sqrt2 = 0.70710678118654752440; /* 1/sqrt(2) */
    static double _sqrtpi = 0.564189583546;        /* 1/sqrt(pi) */

    double x, y, z, k, sum;

    z = fabs(dl - average) / stdv;
    x = z * _sqrt2;
    y = _sqrtpi * exp(-x * x);
    z *= z;
    k = 1.0;
    sum = 0.0;
    do
    {
        sum += x;
        k += 2.0;
        x *= z / k;
    } while (sum + x > sum);
    z = 0.5 - y * sum; /* probability of each tail */
    return (dl > average) ? z : 1.0 / z;
}

void calculate_zscore(double a, int maxdinucleotides, int min, int max, char *filename, int startsite, char *chrom)
{
    static double pideg = 57.29577951; /* 180/pi */
    char *bestantisyn;
    FILE *file;
    unsigned seqlength, i, j;
    int fromdin, todin, din, nucleotides;
    long begintime, endtime;
    double dl, slope, probability, bestdl;
    float initesum;

    fromdin = min;
    todin = max;
    printf("calculating zscore\n");

    file = open_file(1, filename, "");
    if (file == NULL)
    {
        printf("couldn't open %s!\n", filename);
        return;
    }
    seqlength = input_sequence(file, 2 * maxdinucleotides, 0);
    fclose(file);

    file = open_file(0, filename, "Z-SCORE.bedgraph");
    if (file == NULL)
    {
        free(sequence);
        return;
    }

    if (todin > maxdinucleotides)
    {
        todin = maxdinucleotides;
    }
    if (fromdin > todin)
    {
        fromdin = todin;
    }
    nucleotides = 2 * todin;

    a /= 2.0;
    initesum = 10.0 * todin;
    bestantisyn = (char *)malloc(nucleotides + 1);

    time(&begintime);
    for (i = 0; i < seqlength; i++)
    {
        assign_bzenergy_index(nucleotides, sequence + i);
        bestdl = 50.0;
        for (din = fromdin; din <= todin; din++)
        {
            best_esum = initesum;
            deltatwist = a * (double)din;
            antisyn[2 * din] = 0;
            anti_syn_energy(0, din, 0.0); /* esum = 0.0 */
            dl = find_delta_linking(din);
            if (dl < bestdl)
            {
                bestdl = dl;
                strcpy(bestantisyn, best_antisyn);
            }
        }
        probability = assign_probability(bestdl);
        fprintf(file, "%s\t%u\t%u\t%f\n", chrom, startsite + i, startsite + i + 1, probability);
    }
    time(&endtime);
    free(bestantisyn);
    fclose(file);
    printf("\n run time=%ld sec\n", endtime - begintime);
    free(sequence);
}

FILE *open_file(int mode, char *filename, char *typestr)
{
    static char *iostr[] = {"output", "input"};
    static char *rwstr[] = {"w", "r"};
    char *fullfile;

    FILE *file;
    file = NULL;
    fullfile = (char *)malloc(sizeof(char) * (strlen(filename) + strlen(typestr) + 1));
    strcpy(fullfile, filename);
    if (strlen(typestr) != 0)
    {
        strcat(fullfile, ".");
        strcat(fullfile, typestr);
    }
    printf("opening %s\n", fullfile);

    file = fopen(fullfile, rwstr[mode]);
    free(fullfile);
    return file;
}

unsigned input_sequence(FILE *file, int nucleotides, int showfile)
{
    unsigned length, i, j;
    char c;

    printf("inputting sequence\n");

    length = 0; /* count how many bases */
    while (j = 0, fgets(tempstr, 128, file) != NULL)
    {
        while ((c = tempstr[j++]) != 0)
        {
            if (c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'n' || c == 'A' || c == 'T' || c == 'G' || c == 'C' || c == 'N')
            {
                length++;
            }
        }
    }
    sequence = (char *)malloc(length + nucleotides);

    rewind(file);

    if (showfile)
    {
        printf("\n");
    }
    i = 0;
    while (j = 0, fgets(tempstr, 128, file) != NULL)
    {
        while ((c = tempstr[j++]) != 0)
        {
            // Check if c is a nucleotide or 'N' (case insensitive)
            if (strchr("atgcnATGCN", c))
            {
                c = tolower(c);
                sequence[i++] = c;
            }
            else if (strchr("bdefhijklmopqrsuvwxyzBDEFHIJKLMOPQRUVWXYZ0123456789>#$", c))
            {
                // Print error message and exit if c is invalid
                fprintf(stderr, "Error: Invalid character '%c' in input string. Character number = '%u'\nAllowed characters = atgcnATGCN. Exiting. \n", c, i);
                exit(EXIT_FAILURE); // Exit on error
            }
            else
            {
                continue;
            }
        }
    }

    for (j = 0; j < nucleotides; j++) /* assume circular nucleotides */
    {
        sequence[i++] = sequence[j];
    }

    return length;
}

int main(int argc, char *argv[])
{
    static double rt = 0.59004;       /* 0.00198*298 */
    static double a = 0.357, b = 0.4; /* a = 2 * (1/10.5 + 1/12) */
    double ab;
    int i, j, nucleotides, dinucleotides, select;
    int min, max;
    int startsite;
    char *chrom = NULL;
    char *fn = NULL;
    char *extensions[] = {".fasta",".txt", ".fna", ".fa"};
    int num_extensions = sizeof(extensions) / sizeof(extensions[0]);
    char *exch = NULL;
     if (argc == 7)
    {
        chrom = malloc(strlen(argv[5]) + 1);
        strcpy(chrom, argv[5]);
        startsite = atoi((char *)argv[6]);
    }
    else
    {
        printf("usage: z-nome_hunt windowsize minsize maxsize datafile (chromosome_name) (startsite)\n");
        exit(1);
    }

    tempstr = (char *)malloc(128);
    dinucleotides = atoi((char *)argv[1]);
    min = atoi((char *)argv[2]);
    max = atoi((char *)argv[3]);

    printf("dinucleotides %d\n", dinucleotides);
    printf("min/max %d %d\n", min, max);
    printf("chromosome name %s\n", chrom);
    printf("start site %d\n", startsite);

    nucleotides = 2 * dinucleotides;

    antisyn = (char *)malloc(nucleotides + 1);
    best_antisyn = (char *)malloc(nucleotides + 1);

    bzindex = (int *)calloc(dinucleotides, sizeof(int));
    bzenergy = (double *)calloc(dinucleotides, sizeof(double));
    best_bzenergy = (double *)calloc(dinucleotides, sizeof(double));

    bztwist = (double *)calloc(dinucleotides, sizeof(double));

    ab = b + b;
    for (i = 0; i < dinucleotides; i++)
    {
        ab += a;
        bztwist[i] = ab;
    }

    for (i = 0; i < 4; i++)
        for (j = 0; j < 25; j++)
            expdbzed[i][j] = exp(-dbzed[i][j] / rt);

    logcoef = (double *)calloc(dinucleotides, sizeof(double));
    exponent = (double *)calloc(dinucleotides, sizeof(double));
    calculate_zscore(a, dinucleotides, min, max, (char *)argv[4], startsite, (char *)chrom);

    free(exponent);
    free(logcoef);
    free(bztwist);
    free(best_bzenergy);
    free(bzenergy);
    free(bzindex);
    free(best_antisyn);
    free(antisyn);
    free(tempstr);
    return 0;
}
