#include <stdio.h>
#include "dna_analiz.h"
int main(int argc, char *argv[])

{
    if (argc != 3)
    {
        printf("Kullanim: %s <doku_analiz.fasta> <rapor_sonuc.txt>\n", argv[0]);
        return 1;
    }

    char *girdi_dosyasi = argv[1];
    char *cikti_dosyasi = argv[2];

    // Dosyayı analiz etmeden önce ya da içinde gerçek uzunluğu yakalayıp basalım:
    // dosyadanAnalizEt fonksiyonu zaten içeride uzunluğu buluyor.
    // İstersen dosyadanAnalizEt fonksiyonunun içine printf("%d\n", uzunluk); ekleyebilirsin.

    int durum = dosyadanAnalizEt(girdi_dosyasi, cikti_dosyasi);

    if (durum == 0)
    {
        printf("Analiz basariyla tamamlandi. Rapor '%s' dosyasina yazildi.\n", cikti_dosyasi);
    }
    else
    {
        printf("Analiz sirasinda hata olustu!\n");
    }

    return durum;
}