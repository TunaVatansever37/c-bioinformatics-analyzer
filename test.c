#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "dna_sentez/dna_analiz.h"
int main()
{
    int a = 0, t = 0, g = 0, c = 0;

    char yazdir[100] = "";
    char dna[100] = "";
    char rna[100] = "";
    char comp[100] = "";

    printf("======================================\n");
    printf("  MINI BIYOINFOMATIK ANALIZ PROGRAMI  \n");
    printf("======================================\n");
    printf("LUTFEN BIR DNA DIZILIMI GIRINIZ:\n ");

    if (fgets(dna, sizeof(dna), stdin) == NULL)
    {
        printf("HATA: Girdi alinamadi!\n");
        return 1;
    }

    dna[strcspn(dna, "\n")] = '\0';
    int uzunluk = strlen(dna);

    for (int i = 0; i < uzunluk; i++)
    {
        dna[i] = toupper(dna[i]);

        if (dna[i] != 'A' && dna[i] != 'G' && dna[i] != 'C' && dna[i] != 'T')
        {
            printf("HATALI DNA DIZILIMI GIRDINIZ\n");
            return 1;
        }
    }

    dnaDonustur(dna, rna, uzunluk);
    dna_tamamla(dna, comp, uzunluk);
    nukleotid_say(dna, uzunluk, &a, &t, &c, &g);
    float oran = gc_orani(g, c, uzunluk);
    int sonuc = protein_sentezi(rna, uzunluk, yazdir);

    if (sonuc == 0)
    {
        printf("HATA: GECERLI RNA KODONU BULUNAMADI\n");
    }

    printf("\n----------------------\n");
    printf("ANALIZ VERILERI\n");
    printf("----------------------\n");
    printf("ORIJINAL DNA: %s\n", dna);
    printf("RNA SEKANSI: %s\n", rna);
    printf("TERS - TAMAMLAYICI DNA: %s\n", comp);

    printf("\n--- Istatistikler ---\n");
    printf("Toplam nukleotid: %d\n", uzunluk);
    printf("GC Orani: %% %.2f\n", oran);

    printf("\n--- Protein Sentezi ---\n");
    printf("Sentezlenen proteinler:\n");
    printf("%s\n", yazdir);

    return 0;
}