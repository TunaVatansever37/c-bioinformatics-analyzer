#include <stdio.h>
#include <string.h>
#include "dna_analiz.h"

/**
 *Fonksiyon,  dna diziliminde yer alan a(adenin) t(timin) g(guanin) ve c(sitozin) nükleotidlerinin tamamlayıcı nükleotidleri ile eşler.
 *nükleotidlerin tamamlayıcı nükleotidleri: a - t, t - a, g - c, c - g
 * @param dna -  nükleotid dizilimin bulunduğu sarmal yapıdaki yönetici molekül
 * @param comp - nükleotidlerin tamamlayıcı nükleotid karşılıkları
 * @param uzunluk - dna dizilimindeki karakterlerin toplam uzunluğu
 */
void dna_tamamla(char dna[], char comp[], int uzunluk)
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
 *
 * Fonksiyon, verilen DNA karakter dizisini baştan sona tarar ve her bir
 * nükleotitin (A, T, C, G) kaç kez geçtiğini sayarak, adresleri iletilen
 * işaretçi (pointer) değişkenlerine yazar.
 *
 * @param dna     Analiz edilecek DNA karakter dizisi (string)
 * @param uzunluk DNA dizisinin karakter sayısı
 * @param a       'A' (Adenin) sayısının saklanacağı işaretçi
 * @param t       'T' (Timin) sayısının saklanacağı işaretçi
 * @param c       'C' (Sitozin) sayısının saklanacağı işaretçi
 * @param g       'G' (Guanin) sayısının saklanacağı işaretçi
 */
void nukleotid_say(char dna[], int uzunluk, int *a, int *t, int *c, int *g)
{
    /* Dizideki her bir karakteri sırayla kontrol et */
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
            /* Geçersiz karakterler 'gecerliMi' fonksiyonunda süzüldüğü için
               burada ek bir işlem yapmaya gerek yoktur. */
            break;
        }
    }
}

/*
*Fonksiyon, dna düzlemi içerisinde yer alan g(guanin) ve c(sitozin)' in miktarının yzüdesini gösterir
*örneğin: dna dizilimi = (gcta) gc orani = %50
* @param g - dna diziliminde bulunan Guanin nükleotidi
* @param c -  dna diziliminde bulunan sitozin nükleotidi
* @param uzunluk - dna dizilimindeki karakterlerin toplam uzunluğu


*/
float gc_orani(int g, int c, int uzunluk)
{
    if (uzunluk == 0)
        return 0.0f;
    return ((float)(g + c) / (float)uzunluk) * 100.0f;
}

/**
 * Fonksiyon, dna düzlemi içerisinde yer alan t(Timin) nükleotidini rna diziliminine uygun yapmak için u(urasil) nükleotidine çevirir.
 * örneğin: dna dizilimi = tagcgtacg rna dizilimli yeni hali = uagcguacg
 * bu sayede dna diziliminin rna dizilimine çevirdiğimizde nasıl göründüğünü öğreniriz.
 * @param dna - nükleotid dizilimin bulunduğu sarmal yapıdaki yönetici molekül
 * @param rna - dna dizilimindeki timin nükleotidlerinin urasile çevrilmiş hali
 * @param uzunluk -  dna dizilimindeki karakterlerin toplam uzunluğu
 */
void dnaDonustur(char *dna, char *rna, int uzunluk)
{
    for (int i = 0; i < uzunluk; i++)
    {
        rna[i] = (dna[i] == 'T') ? 'U' : dna[i];
    }
    rna[uzunluk] = '\0';
}

/**
 * Fonksiyon, rna dizilimi içerisinde protein sentezi yapabilecek 3 er li nükleotid dizilimlerini alıp hangi sentezi yapabileceklerini söyler
 * örneğin - gcutata - gcu üçlüsü Ala(Alanin) sentezi yapabilir
 * @param rna - dna dizilimindeki timin nükleotidlerinin urasile çevrilmiş hali
 * @param uzunluk - dna dizilimindeki karakterlerin toplam uzunluğu
 * @param yazdir - strcat ile üzerine kopya yazacağımız boş alan
 */
int protein_sentezi(char rna[], int uzunluk, char yazdir[])
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
            strcat(yazdir, "-Asp");
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
        else if ((rna[i] == 'U' && rna[i + 1] == 'A' && rna[i + 2] == 'Y') || (rna[i] == 'U' && rna[i + 1] == 'A' && rna[i + 2] == 'C'))
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