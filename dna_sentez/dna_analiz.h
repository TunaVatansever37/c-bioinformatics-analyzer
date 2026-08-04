#ifndef DNA_ANALIZ_H
#define DNA_ANALIZ_H

int protein_sentezi(char rna[], int uzunluk, char yazdir[]);
float gc_orani(int g, int c, int uzunluk);
void nukleotid_say(char dna[], int uzunluk, int *a, int *t, int *c, int *g);
void dna_tamamla(char dna[], char comp[], int uzunluk);
void dnaDonustur(char *dna, char *rna, int uzunluk);

#endif