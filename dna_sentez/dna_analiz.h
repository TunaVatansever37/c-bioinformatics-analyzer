#ifndef DNA_ANALIZ_H
#define DNA_ANALIZ_H

int proteinSentezi(char rna[], int uzunluk, char yazdir[]);
float gcOrani(int g, int c, int uzunluk);
void nukleotidSay(char dna[], int uzunluk, int *a, int *t, int *c, int *g);
void dnaTamamla(char dna[], char comp[], int uzunluk);
void dnaDonustur(char *dna, char *rna, int uzunluk);
int dosyadanAnalizEt(char *input, char *output);

#endif