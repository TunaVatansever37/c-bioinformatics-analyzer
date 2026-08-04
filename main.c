#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "dna_sentez/dna_analiz.h"

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Kullanim: ./bioinfo_app <girdi_dosyasi.fasta> <cikti_raporu.txt>\n");
        return 1;
    }

    char *girdi_dosyasi = argv[1];
    char *cikti_dosyasi = argv[2];

       FILE *fp_in = fopen(girdi_dosyasi, "r");
    if (fp_in == NULL)
    {
        printf("Hata: Girdi dosyasi acilamadi!\n");
        return 1;
    }

    char satir[256];
    char baslik[256] = "";
    char dna[512] = "";

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
        }
    }
    fclose(fp_in);

    int uzunluk = strlen(dna);
    if (uzunluk == 0)
    {
        printf("Hata: Dosyada DNA dizilimi bulunamadi!\n");
        return 1;
    }

    for (int i = 0; i < uzunluk; i++)
    {
        dna[i] = toupper(dna[i]);
        if (dna[i] != 'A' && dna[i] != 'G' && dna[i] != 'C' && dna[i] != 'T')
        {
            printf("HATALI DNA DIZILIMI!\n");
            return 1;
        }
    }

    int a = 0, t = 0, g = 0, c = 0;
    char yazdir[256] = "";
    char rna[512] = "";
    char comp[512] = "";

    dnaDonustur(dna, rna, uzunluk);
    dnaTamamla(dna, comp, uzunluk);
    nukleotidSay(dna, uzunluk, &a, &t, &c, &g);
    float oran = gcOrani(g, c, uzunluk);
    proteinSentezi(rna, uzunluk, yazdir);

    FILE *fp_out = fopen(cikti_dosyasi, "w");
    if (fp_out == NULL)
    {
        printf("Hata: Rapor dosyasi olusturulamadi!\n");
        return 1;
    }

    fprintf(fp_out, "==================================================\n");
    fprintf(fp_out, "BIYOINFORMATIK ANALIZ RAPORU\n");
    fprintf(fp_out, "==================================================\n");
    fprintf(fp_out, "Sekans Basligi : %s\n", baslik);
    fprintf(fp_out, "------------------ DIZILIMLER --------------------\n");
    fprintf(fp_out, "Orijinal DNA (5'->3') : %s\n", dna);
    fprintf(fp_out, "RNA Sekansi (5'->3') : %s\n", rna);
    fprintf(fp_out, "Ters-Tamamlayici DNA : %s\n", comp);
    fprintf(fp_out, "----------------- ISTATISTIKLER ------------------\n");
    fprintf(fp_out, "Toplam Nukleotid : %d\n", uzunluk);
    fprintf(fp_out, "A Sayisi : %d\n", a);
    fprintf(fp_out, "T Sayisi : %d\n", t);
    fprintf(fp_out, "G Sayisi : %d\n", g);
    fprintf(fp_out, "C Sayisi : %d\n", c);
    fprintf(fp_out, "GC Orani : %% %.2f\n", oran);
    fprintf(fp_out, "------------- PROTEIN TRANSLASYONU ---------------\n");
    fprintf(fp_out, "Sentezlenen Protein : %s\n", yazdir);
    fprintf(fp_out, "==================================================\n");

    fclose(fp_out);

    printf("Analiz basariyla tamamlandi. Rapor '%s' dosyasina yazildi.\n", cikti_dosyasi);
    return 0;
}