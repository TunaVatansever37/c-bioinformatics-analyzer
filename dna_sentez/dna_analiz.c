#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "dna_analiz.h"

int dosyadanAnalizEt(char *input, char *output)
{
    FILE *fp_in = fopen(input, "r");
    if (fp_in == NULL)
    {
        printf("HATA: Girdi dosyasi (%s) acilamadi\n", input);
        return 1;
    }

    fseek(fp_in, 0, SEEK_END);
    long dosya_boyutu = ftell(fp_in);
    rewind(fp_in);

    if (dosya_boyutu <= 0)
    {
        printf("Hata: Dosya bos!\n");
        fclose(fp_in);
        return 1;
    }

    char *dna = malloc((dosya_boyutu + 1) * sizeof(char));
    if (dna == NULL)
    {
        printf("Hata: Bellek tahsis edilemedi!\n");
        fclose(fp_in);
        return 1;
    }
    dna[0] = '\0';

    char satir[256];
    char baslik[256] = "";
    size_t dna_uzunlugu = 0;

    while (fgets(satir, sizeof(satir), fp_in) != NULL)
    {
        satir[strcspn(satir, "\r\n")] = 0;
        if (strlen(satir) == 0)
            continue;

        if (satir[0] == '>')
        {
            strcpy(baslik, satir + 1);
        }
        else
        {
            strcat(dna, satir);
            dna_uzunlugu += strlen(satir);
        }
    }
    fclose(fp_in);

    int uzunluk = strlen(dna);
    if (uzunluk == 0)
    {
        printf("Hata: Dosyada DNA dizilimi bulunamadi!\n");
        free(dna);
        return 1;
    }

    for (int i = 0; i < uzunluk; i++)
    {
        dna[i] = toupper(dna[i]);
    }

    int a = 0, t = 0, g = 0, c = 0;
    char *rna = malloc((uzunluk + 1) * sizeof(char));
    char *comp = malloc((uzunluk + 1) * sizeof(char));
    char *yazdir = malloc((uzunluk + 1) * sizeof(char));

    if (rna == NULL || comp == NULL || yazdir == NULL)
    {
        printf("Hata: Ek bellek alanlari olusturulamadi!\n");
        free(dna);
        free(rna);
        free(comp);
        free(yazdir);
        return 1;
    }

    dnaDonustur(dna, rna, uzunluk);
    dnaTamamla(dna, comp, uzunluk);
    nukleotidSay(dna, uzunluk, &a, &t, &c, &g);
    float oran = gcOrani(g, c, uzunluk);
    proteinSentezi(rna, uzunluk, yazdir);

    FILE *fp_out = fopen(output, "w");
    if (fp_out == NULL)
    {
        printf("HATA: cikti dosyasi (%s) olusturulamadi!\n", output);
        free(dna);
        free(rna);
        free(comp);
        free(yazdir);
        return 1;
    }

    fprintf(fp_out, "==================================================\n");
    fprintf(fp_out, "                BIYOINFORMATIK ANALIZ RAPORU\n");
    fprintf(fp_out, "==================================================\n");
    fprintf(fp_out, "Sekans Basligi        : %s\n", baslik);
    fprintf(fp_out, "------------------ DIZILIMLER --------------------\n");
    fprintf(fp_out, "Orijinal DNA (5'->3') : %s\n", dna);
    fprintf(fp_out, "RNA Sekansi  (5'->3') : %s\n", rna);
    fprintf(fp_out, "Ters-Tamamlayici DNA  : %s\n", comp);
    fprintf(fp_out, "----------------- ISTATISTIKLER ------------------\n");
    fprintf(fp_out, "Toplam Nukleotid      : %d\n", uzunluk);
    fprintf(fp_out, "A Sayisi              : %d\n", a);
    fprintf(fp_out, "T Sayisi              : %d\n", t);
    fprintf(fp_out, "G Sayisi              : %d\n", g);
    fprintf(fp_out, "C Sayisi              : %d\n", c);
    fprintf(fp_out, "GC Orani              : %% %.2f\n", oran);
    fprintf(fp_out, "------------- PROTEIN TRANSLASYONU ---------------\n");
    fprintf(fp_out, "Sentezlenen Protein   : %s\n", yazdir);
    fprintf(fp_out, "==================================================\n");

    fclose(fp_out);
    free(dna);
    free(rna);
    free(comp);
    free(yazdir);

    return 0;
}

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
            strcat(yazdir, "Met-");
            bulundu = 1;
        }
        else if ((rna[i] == 'G' && rna[i + 1] == 'C' && rna[i + 2] == 'U') || (rna[i] == 'G' && rna[i + 1] == 'C' && rna[i + 2] == 'C') || (rna[i] == 'G' && rna[i + 1] == 'C' && rna[i + 2] == 'A') || (rna[i] == 'G' && rna[i + 1] == 'C' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Ala-");
            bulundu = 1;
        }
        else if ((rna[i] == 'C' && rna[i + 1] == 'G' && rna[i + 2] == 'U') || (rna[i] == 'C' && rna[i + 1] == 'G' && rna[i + 2] == 'C') || (rna[i] == 'C' && rna[i + 1] == 'G' && rna[i + 2] == 'A') || (rna[i] == 'C' && rna[i + 1] == 'G' && rna[i + 2] == 'G') || (rna[i] == 'A' && rna[i + 1] == 'G' && rna[i + 2] == 'A') || (rna[i] == 'A' && rna[i + 1] == 'G' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Arg-");
            bulundu = 1;
        }
        else if ((rna[i] == 'A' && rna[i + 1] == 'A' && rna[i + 2] == 'U') || (rna[i] == 'A' && rna[i + 1] == 'A' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Asn-");
            bulundu = 1;
        }
        else if ((rna[i] == 'G' && rna[i + 1] == 'A' && rna[i + 2] == 'U') || (rna[i] == 'G' && rna[i + 1] == 'A' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Asp-");
            bulundu = 1;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'A' && rna[i + 2] == 'A') || (rna[i] == 'U' && rna[i + 1] == 'A' && rna[i + 2] == 'G') || (rna[i] == 'U' && rna[i + 1] == 'G' && rna[i + 2] == 'A'))
        {
            strcat(yazdir, "[STOP]");
            break;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'G' && rna[i + 2] == 'U') || (rna[i] == 'U' && rna[i + 1] == 'G' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Cys-");
            bulundu = 1;
        }
        else if ((rna[i] == 'C' && rna[i + 1] == 'A' && rna[i + 2] == 'A') || (rna[i] == 'C' && rna[i + 1] == 'A' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Gln-");
            bulundu = 1;
        }
        else if ((rna[i] == 'G' && rna[i + 1] == 'A' && rna[i + 2] == 'A') || (rna[i] == 'G' && rna[i + 1] == 'A' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Glu-");
            bulundu = 1;
        }
        else if ((rna[i] == 'G' && rna[i + 1] == 'G' && rna[i + 2] == 'U') || (rna[i] == 'G' && rna[i + 1] == 'G' && rna[i + 2] == 'C') || (rna[i] == 'G' && rna[i + 1] == 'G' && rna[i + 2] == 'A') || (rna[i] == 'G' && rna[i + 1] == 'G' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Gly-");
            bulundu = 1;
        }
        else if ((rna[i] == 'C' && rna[i + 1] == 'A' && rna[i + 2] == 'U') || (rna[i] == 'C' && rna[i + 1] == 'A' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "His-");
            bulundu = 1;
        }
        else if ((rna[i] == 'A' && rna[i + 1] == 'U' && rna[i + 2] == 'U') || (rna[i] == 'A' && rna[i + 1] == 'U' && rna[i + 2] == 'C') || (rna[i] == 'A' && rna[i + 1] == 'U' && rna[i + 2] == 'A'))
        {
            strcat(yazdir, "Ile-");
            bulundu = 1;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'G' && rna[i + 2] == 'G'))
        {
            strcat(yazdir, "Trp-");
            bulundu = 1;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'A' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Tyr-");
            bulundu = 1;
        }
        else if ((rna[i] == 'U' && rna[i + 1] == 'U' && rna[i + 2] == 'U') || (rna[i] == 'U' && rna[i + 1] == 'U' && rna[i + 2] == 'C'))
        {
            strcat(yazdir, "Phe-");
            bulundu = 1;
        }
    }
    return bulundu;
}