#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "dna_analiz.h"

/**
 * Fonksiyon, dna diziliminde yer alan a(adenin) t(timin) g(guanin) ve c(sitozin) nükleotidlerinin tamamlayıcı nükleotidleri ile eşler.
 */
void dnaTamamla(char dna[], char comp[], int uzunluk)
{
    for (int i = 0; i < uzunluk; i++)
    {
        switch (dna[i])
        {
        case 'A':
            comp[i] = 'T';
            break;
        case 'T':
            comp[i] = 'A';
            break;
        case 'C':
            comp[i] = 'G';
            break;
        case 'G':
            comp[i] = 'C';
            break;
        }
    }
    comp[uzunluk] = '\0';
}

/**
 * @brief DNA dizisindeki A, T, C ve G nükleotitlerinin sayılarını hesaplar.
 */
void nukleotidSay(char dna[], int uzunluk, int *a, int *t, int *c, int *g)
{
    for (int i = 0; i < uzunluk; i++)
    {
        switch (dna[i])
        {
        case 'A':
            (*a)++;
            break;
        case 'C':
            (*c)++;
            break;
        case 'G':
            (*g)++;
            break;
        case 'T':
            (*t)++;
            break;
        default:
            break;
        }
    }
}

float gcOrani(int g, int c, int uzunluk)
{
    if (uzunluk == 0)
        return 0.0f;
    return ((float)(g + c) / (float)uzunluk) * 100.0f;
}

void dnaDonustur(char *dna, char *rna, int uzunluk)
{
    for (int i = 0; i < uzunluk; i++)
    {
        rna[i] = (dna[i] == 'T') ? 'U' : dna[i];
    }
    rna[uzunluk] = '\0';
}

int proteinSentezi(char rna[], int uzunluk, char yazdir[])
{
    yazdir[0] = '\0';
    int bulundu = 0;

    for (int i = 0; i < uzunluk - 2; i += 3)
    {
        if (rna[i] == 'A' && rna[i + 1] == 'U' && rna[i + 2] == 'G')
        {
            strcat(yazdir, "Met");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'G' && rna[i + 1] == 'C' && rna[i + 2] == 'U') || (rna[i] == 'G' && rna[i + 1] == 'C' && rna[i + 2] == 'C') || (rna[i] == 'G' && rna[i + 1] == 'C' && rna[i + 2] == 'A') || (rna[i] == 'G' && rna[i + 1] == 'C' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Ala");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'C' && rna[i + 1] == 'G' && rna[i + 2] == 'U') || (rna[i] == 'C' && rna[i + 1] == 'G' && rna[i + 2] == 'C') || (rna[i] == 'C' && rna[i + 1] == 'G' && rna[i + 2] == 'A') || (rna[i] == 'C' && rna[i + 1] == 'G' && rna[i + 2] == 'G') || (rna[i] == 'A' && rna[i + 1] == 'G' && rna[i + 2] == 'A') || (rna[i] == 'A' && rna[i + 1] == 'G' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Arg");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'A' && rna[i + 1] == 'A' && rna[i + 2] == 'U') || (rna[i] == 'A' && rna[i + 1] == 'A' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Asn");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'G' && rna[i + 1] == 'A' && rna[i + 2] == 'U') || (rna[i] == 'G' && rna[i + 1] == 'A' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Asp");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'A' && rna[i + 2] == 'A') || (rna[i] == 'U' && rna[i + 1] == 'A' && rna[i + 2] == 'G') || (rna[i] == 'U' && rna[i + 1] == 'G' && rna[i + 2] == 'A'))
        {
            strcat(yazdir, "[STOP]");
            break;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'G' && rna[i + 2] == 'U') || (rna[i] == 'U' && rna[i + 1] == 'G' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Cys");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'C' && rna[i + 1] == 'A' && rna[i + 2] == 'A') || (rna[i] == 'C' && rna[i + 1] == 'A' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Gln");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'G' && rna[i + 1] == 'A' && rna[i + 2] == 'A') || (rna[i] == 'G' && rna[i + 1] == 'A' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Glu");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'G' && rna[i + 1] == 'G' && rna[i + 2] == 'U') || (rna[i] == 'G' && rna[i + 1] == 'G' && rna[i + 2] == 'C') || (rna[i] == 'G' && rna[i + 1] == 'G' && rna[i + 2] == 'A') || (rna[i] == 'G' && rna[i + 1] == 'G' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Gly");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'C' && rna[i + 1] == 'A' && rna[i + 2] == 'U') || (rna[i] == 'C' && rna[i + 1] == 'A' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "His");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'A' && rna[i + 1] == 'U' && rna[i + 2] == 'U') || (rna[i] == 'A' && rna[i + 1] == 'U' && rna[i + 2] == 'C') || (rna[i] == 'A' && rna[i + 1] == 'U' && rna[i + 2] == 'A'))
        {
            strcat(yazdir, "Ile");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'G' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Trp");
            strcat(yazdir, "-");
            bulundu = 1;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'A' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Tyr");
            bulundu = 1;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'U' && rna[i + 2] == 'U') || (rna[i] == 'U' && rna[i + 1] == 'U' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Phe");
            strcat(yazdir, "-");
            bulundu = 1;
        }
    }
    return bulundu;
}